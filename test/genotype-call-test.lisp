;;;
;;; Copyright (c) 2011-2012 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of genotype-call.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.genotype-call-test)

(defparameter *possible-alleles*
  '("AA" "AC" "AG" "AT" "CC" "CG" "GC" "GG" "TA" "TC" "TG" "TT" "DI" "ID" "NA"))

(defparameter *possible-strands*
  '("Bot" "BOT" "M" "MINUS" "P" "PLUS" "Top" "TOP"))

(defparameter *example-snp-counts*
    '(("16" 4) ("2" 4) ("9" 2)))

;; To check the SIM file contents:
;;
;;  hexdump -v -s 16 -e '1/255 "%s" ": " 20/4 " %0.12f" "\n"' example.sim
;;
(defparameter *sim-intensities-hexdump*
  '(0.034164227545 0.074980102479
                   0.180595964193 0.046113856137
                   1.067971110344 0.008573265746
                   0.891155481339 0.873829364777
                   0.600366592407 0.002887635026
                   1.005514144897 0.987805783749
                   -0.000895529287 1.537799477577
                   0.014101119712 0.762843370438
                   0.712541103363 0.023448554799
                   0.901769757271 1.160114407539))

(defun round-float (f)
  (let ((fact (expt 10 9)))
    (/ (fround (* fact f)) fact)))

(defun compare-intensities (expected observed)
  (ensure (= (length expected) (length observed))
          :report "intensity sequences were not the same length: ~a ~a"
          :arguments (expected observed))
  (mapcar (lambda (x y)
            (let ((rx (round-float x))
                  (ry (round-float y)))
              (ensure (= rx ry)
                      :report "intensity comparison failed: ~a (rounded to ~a) != ~a (rounded to ~a)"
                      :arguments (x rx y ry))))
          expected (coerce observed 'list)))


(deftestsuite genotype-call-tests ()
  ())

(addtest (genotype-call-tests) read-bpm/1
  (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
    (let ((manifest (read-bpm stream))
          (chrs (mapcar #'first *example-snp-counts*))
          (counts (mapcar #'second *example-snp-counts*)))
      (ensure (= 10 (num-snps-of manifest)))
      (ensure (equal chrs (chromosomes-of manifest)))
      (ensure (equal counts (mapcar (lambda (chr)
                                      (num-snps-of manifest
                                                   :key #'snp-chromosome
                                                   :test (lambda (x)
                                                           (string= chr x))))
                                    chrs))))))

(addtest (genotype-call-tests) read-bpm/2
  (with-open-file (stream (merge-pathnames "data/example_unsorted.bpm.csv"))
    (ensure (read-bpm stream :strict-ordering nil)))
  (with-open-file (stream (merge-pathnames "data/example_unsorted.bpm.csv"))
    (ensure-condition (malformed-file-error)
      (read-bpm stream :strict-ordering t))))

(addtest (genotype-call-tests) normalize-alleles/1
  (dolist (alleles *possible-alleles*)
    (dolist (strand (remove-duplicates
                     (mapcar (lambda (str)
                               (genotype-call::parse-strand nil str))
                             *possible-strands*)))
      (let ((norm (genotype-call::normalize-alleles alleles strand)))
        (cond ((string= "NA" norm)
               nil)
              ((eql #\B strand)
               (ensure (string= (map 'simple-string
                                     #'genotype-call::complement-allele
                                     alleles) norm)))
              (t
               (string= alleles norm)))))))

(addtest (genotype-call-tests) valid-alleles-p/1
  (ensure (every #'valid-alleles-p *possible-alleles*)))

(addtest (genotype-call-tests) make-chromosome-p/1
  (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
    (let ((manifest (read-bpm stream))
          (chrs (mapcar #'first *example-snp-counts*)))
      (dolist (chr chrs)
        (let ((pred (make-chromosome-p manifest chr #'string=)))
          (ensure (funcall pred chr))
          (ensure (not (funcall pred "NOT A VALID CHROMOSOME"))))))))

(addtest (genotype-call-tests) gtc-open/close/1
  (with-open-file (stream (merge-pathnames "data/example.gtc")
                          :element-type 'octet)
    (let ((gtc (gtc-open stream)))
      (ensure gtc)
      (ensure (= 3 (version-of gtc)))
      (ensure (= 24 (length (toc-of gtc))))
      (ensure (gtc-close gtc)))))

(addtest (genotype-call-tests) with-gtc/1
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
     (ensure (= 3 (version-of gtc)))
     (ensure (= 24 (length (toc-of gtc))))))

(addtest (genotype-call-tests) data-field-of/1
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
    (ensure (= 660447 (data-field-of gtc :num-snps)))
    (ensure (string= "test_sample" (data-field-of gtc :sample-name)))
    (ensure (string= "test_plate_01" (data-field-of gtc :sample-plate)))
    (ensure (string= "A01" (data-field-of gtc :sample-well)))
    (ensure (string= "Human670-QuadCustom_v1_A.egt"
                     (data-field-of gtc :cluster-file)))
    (ensure (string= "Human670-QuadCustom_v1_A.bpm"
                     (data-field-of gtc :snp-manifest)))
    (ensure (string= "Tuesday, June 14, 2011 1:55:35 PM"
                     (data-field-of gtc :imaging-date)))
    (ensure (string= "6/14/2011 2:12 PM"
                     (data-field-of gtc :autocall-date)))
    (ensure (string= "1.6.2.2"
                     (data-field-of gtc :autocall-version)))

    (ensure (= 24 (length (data-field-of gtc :normalization-xforms))))

    (ensure-condition (invalid-argument-error)
      (data-field-of gtc :invalid-field))
      
    (ensure (= 264 (length (data-field-of gtc :x-controls))))
    (ensure (= 264 (length (data-field-of gtc :y-controls))))
      
    (ensure (= 660447 (length (data-field-of gtc :x-intensities))))
    (ensure (= 660447 (length (data-field-of gtc :y-intensities))))
    (ensure (= 660447 (length (data-field-of gtc :genotypes))))
    (ensure (= 660447 (length (data-field-of gtc :basecalls))))
    (ensure (= 660447 (length (data-field-of gtc :genotype-scores))))))

(addtest (genotype-call-tests) data-field-of/2
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
    (let ((xforms (data-field-of gtc :normalization-xforms)))
      (loop
         for xform across xforms
         do (ensure (every #'consp (mapcar (lambda (field)
                                             (assoc field xform))
                                           *xform-fields*)))))))

(addtest (genotype-call-tests) data-field-of/3
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
    (ensure (every (lambda (x)
                     (and (vectorp x) (integerp (aref x 0))))
                   (mapcar (lambda (field)
                             (data-field-of gtc field))
                           '(:x-controls :y-controls
                             :x-intensities :y-intensities))))))

(addtest (genotype-call-tests) data-field-of/4
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
    (ensure (every (lambda (x)
                     (and (vectorp x) (stringp (aref x 0))))
                   (mapcar (lambda (field)
                             (data-field-of gtc field))
                           '(:genotypes :basecalls))))))

(addtest (genotype-call-tests) data-field-of/5
  (with-gtc (gtc (merge-pathnames "data/example.gtc"))
    (let ((scores (data-field-of gtc :genotype-scores)))
      (ensure (and (vectorp scores) (floatp (aref scores 0)))))))

(addtest (genotype-call-tests) gtc-to-sim/1
  (handler-bind ((error #'leave-tmp-pathname))
    (with-tmp-pathname (tmp :tmpdir (merge-pathnames "data") :type "sim")
      (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
        (let ((manifest (read-bpm stream))
              (count 5))
          (with-sim (sim tmp :direction :output :if-exists :supersede
                         :if-does-not-exist :create)
            (dotimes (n count)
              (with-gtc (gtc (merge-pathnames "data/example.gtc"))
                (ensure (copy-intensities gtc sim manifest)))
              (ensure (= (1+ n) (num-samples-of sim))
                      :report "expected ~d, but found ~d"
                      :arguments ((1+ n) (num-samples-of sim)))))
          (with-sim (sim tmp)
            (ensure (= count (num-samples-of sim)))))))))

(addtest (genotype-call-tests) sim-open-close/1
  (with-open-file (stream (merge-pathnames "data/example.sim")
                          :element-type 'octet)
    (let ((sim (sim-open stream)))
      (ensure sim)
      (ensure (= 1 (version-of sim)))
      (ensure (sim-close sim)))))

(addtest (genotype-call-tests) sim-closed-stream/1
  (ensure-condition (invalid-operation-error)
    (make-instance 'sim :stream (with-open-file
                                    (stream (merge-pathnames "data/example.sim")
                                            :element-type 'octet)
                                  stream))))

(addtest (genotype-call-tests) with-sim/1
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (ensure sim)
    (ensure (= 1 (version-of sim)))
    (ensure (= 255 (name-size-of sim)))
    (ensure (= 5 (num-samples-of sim)))
    (ensure (= 10 (num-probes-of sim)))
    (ensure (= 2 (num-channels-of sim)))
    (ensure (eql 'single-float (format-of sim)))))

(addtest (genotype-call-tests) read-intensities/1
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (dotimes (n (num-samples-of sim))
      (multiple-value-bind (intensities sample-name)
          (read-intensities sim)
        (compare-intensities *sim-intensities-hexdump* intensities)
        (ensure (string= "test_sample" sample-name))))
    ;; Should be at eof
    (ensure (null (read-byte (stream-of sim) nil nil)))))

(addtest (genotype-call-tests) read-intensities/2
  ;; :start argument only
  (loop
     for i from 0 below 10
     do (with-sim (sim (merge-pathnames "data/example.sim"))
          (compare-intensities (subseq *sim-intensities-hexdump* (* 2 i))
                               (read-intensities sim :start i))))
  ;; :end argument only
  (loop
     for i from 0 below 10
     do (with-sim (sim (merge-pathnames "data/example.sim"))
          (compare-intensities (subseq *sim-intensities-hexdump* 0 (* 2 i))
                               (read-intensities sim :end i))))
  ;; both :start and :end arguments
  (loop
     for i from 0 below 8
     for j from 2 below 10
     do (with-sim (sim (merge-pathnames "data/example.sim"))
          (compare-intensities
           (subseq *sim-intensities-hexdump* (* 2 i) (* 2 j))
           (read-intensities sim :start i :end j)))))

(addtest (genotype-call-tests) read-intensities/3
  ;; single intensity channel set
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (ensure (= (num-channels-of sim)
               (length (read-intensities sim :start 0 :end 1)))))
  ;; start < 0
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (ensure-condition (invalid-argument-error)
      (read-intensities sim :start -1)))
  ;; end == start
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (let ((intensities (read-intensities sim :start 0 :end 0)))
      (ensure (vectorp intensities))
      (ensure (zerop (length intensities)))))
   ;; end < start
   (with-sim (sim (merge-pathnames "data/example.sim"))
     (ensure-condition (invalid-argument-error)
       (read-intensities sim :start 0 :end -1))))

(addtest (genotype-call-tests) sim-to-illuminus/1
  (flet ((read-lines (filespec)
           (with-open-file (stream filespec)
             (loop
                for line = (read-line stream nil nil)
                while line
                collect line))))
    (handler-bind ((error #'leave-tmp-pathname))
      (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
        (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
        (let ((manifest (read-bpm stream))
              (count 5))
          (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                         :if-does-not-exist :create)
            (dotimes (n count)
              (with-gtc (gtc (merge-pathnames "data/example.gtc"))
                (copy-intensities gtc sim manifest))))
          (with-tmp-pathname (tmp-iln :tmpdir (merge-pathnames "data")
                                      :type "iln")
            (let ((iln (sim-to-illuminus tmp-iln manifest tmp-sim)))
              (ensure iln))
            (let ((expected (read-lines (merge-pathnames
                                         "data/example.iln")))
                  (observed (read-lines tmp-iln)))
              (ensure (equalp observed expected)
                      :report "expected ~a but found ~a"
                      :arguments (expected observed))))))))))
