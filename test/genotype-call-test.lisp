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
  '(#\A #\C #\G #\T #\D #\I))

(defparameter *possible-allele-pairs*
  (let ((pairs (list '(#\D #\I) '(#\I #\D))))
    (dolist (a '(#\A #\C #\G #\T))
      (dolist (b '(#\A #\C #\G #\T))
        (push (list a b) pairs)))))

(defparameter *possible-strands*
  '("Bot" "BOT" "M" "MINUS" "P" "PLUS" "Top" "TOP"))

(defparameter *example-bpm*
  (load-bpm (merge-pathnames "data/example.bpm.csv") :name "example.bpm"))

(defparameter *alternative-bpm*
  (load-bpm (merge-pathnames "data/alternative.bpm.csv")
            :name "alternative.bpm"))

(defparameter *example-normalized-sim*
  (merge-pathnames "data/example.normalized.sim"))

(defparameter *example-raw-sim*
  (merge-pathnames "data/example.raw.sim"))


;; To check the SIM file contents:
;;
;;  hexdump -v -s 16 -e '1/255 "%s" ": " 20/4 " %0.12f" "\n"' \
;;    example.normalized.sim
;;
(defparameter *sim-intensities-hexdump*
  '(0.034164227545 0.074980102479
    0.180595964193 0.046113856137
    1.067971110344 0.008573265746
    0.891155481339 0.873829364777
    0.600366592407 0.002887635026
    1.005514144897 0.987805783749
    0.000000000000 1.537799477577
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

(defun binary-file= (x y)
  (with-open-file (s1 x :element-type 'octet)
    (with-open-file (s2 y :element-type 'octet)
      (let ((v1 (make-array (file-length s1) :element-type 'octet))
            (v2 (make-array (file-length s2) :element-type 'octet)))
        (read-sequence v1 s1)
        (read-sequence v2 s2)
        (equalp v1 v2)))))

(defun ensure-lines-equal (expected observed)
  (flet ((read-lines (filespec)
           (with-open-file (stream filespec)
             (loop
                for line = (read-line stream nil nil)
                while line
                collect line))))
    (let ((e (read-lines expected))
          (o (read-lines observed)))
      (ensure (equal e o)
              :report "expected ~a but found ~a"
              :arguments (e o)))))

(defun example-gtc-specs ()
  (mapcar (lambda (file uri gender)
            (pairlis '(:result :uri :gender)
                     (list file (puri:uri uri) gender)))
          '("data/example_0000.gtc"
            "data/example_0001.gtc"
            "data/example_0002.gtc"
            "data/example_0003.gtc"
            "data/example_0004.gtc")
          '("urn:wtsi:example_0000"
            "urn:wtsi:example_0001"
            "urn:wtsi:example_0002"
            "urn:wtsi:example_0003"
            "urn:wtsi:example_0004")
          '(0 1 1 2 2)))

(deftestsuite genotype-call-tests ()
  ())

(addtest (genotype-call-tests) read-bpm/1
  (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
    (let ((manifest (read-bpm stream))
          (chrs (mapcar (lambda (x)
                          (format nil "~d" x)) (iota 10 1))))
      (ensure (= 10 (num-snps-of manifest)))
      (let ((observed (chromosomes-of manifest))
            (expected chrs))
        (ensure (equal expected observed)
                :report "expected ~a but found ~a"
                :arguments (expected observed)))
      (let ((observed (mapcar (lambda (chr)
                                (num-snps-of manifest
                                             :key #'snp-chromosome
                                             :test (lambda (x)
                                                     (string= chr x))))
                              chrs))
            (expected (loop repeat 10 collect 1)))
        (ensure (equal expected observed)
                :report "expected ~a but found ~a"
                :arguments (expected observed))))))

(addtest (genotype-call-tests) read-bpm/2
  (with-open-file (stream (merge-pathnames "data/example_unsorted.bpm.csv"))
    (ensure (read-bpm stream :strict-ordering nil)))
  (with-open-file (stream (merge-pathnames "data/example_unsorted.bpm.csv"))
    (ensure-condition (malformed-file-error)
      (read-bpm stream :strict-ordering t))))

(addtest (genotype-call-tests) location</1
  (ensure (not (location< (make-snp 1 "x" "1" 1000 #\A #\A #\+ #\+ 0)
                          (make-snp 1 "x" "1" 1000 #\A #\A #\+ #\+ 0))))
  (ensure (location< (make-snp 1 "x" "1" 1000 #\A #\A #\+ #\+ 0)
                     (make-snp 1 "x" "1" 1001 #\A #\A #\+ #\+ 0)))
  (ensure (not (location< (make-snp 1 "x" "1" 1001 #\A #\A #\+ #\+ 0)
                          (make-snp 1 "x" "1" 1000 #\A #\A #\+ #\+ 0))))
  (ensure (location< (make-snp 1 "x" "2" 1000 #\A #\A #\+ #\+ 0)
                     (make-snp 1 "x" "10" 1000 #\A #\A #\+ #\+ 0)))
  (ensure (location< (make-snp 1 "x" "2" 1001 #\A #\A #\+ #\+ 0)
                     (make-snp 1 "x" "10" 1000 #\A #\A #\+ #\+ 0))))

(addtest (genotype-call-tests) chromosome-boundaries/1
  (mapc (lambda (chr start end)
          (multiple-value-bind (s e)
              (chromosome-boundaries *example-bpm* chr)
            (ensure (= start s)
                    :report "expected start ~a, but found ~a"
                    :arguments (start s))
            (ensure (= end e)
                    :report "expected end ~a, but found ~a"
                    :arguments (end e))))
        (chromosomes-of *example-bpm*)
        (iota 10 0)
        (iota 10 1)))

(addtest (genotype-call-tests) chromosome-boundaries/2
  (ensure-condition invalid-argument-error
    (chromosome-boundaries *example-bpm* "1" :key #'snp-chromosome
                           :test (lambda (x)
                                   (string= x "no such chromosome")))))

(addtest (genotype-call-tests) normalize-allele/1
  (dolist (allele *possible-alleles*)
    (dolist (strand (remove-duplicates
                     (mapcar (lambda (str)
                               (genotype-call::parse-strand nil str))
                             *possible-strands*)))
      (let ((norm (genotype-call::normalize-allele allele strand)))
        (cond ((char= #\. norm)
               nil)
              ((eql #\B strand)
               (ensure (char= (genotype-call::complement-allele allele) norm)))
              (t
               (char= allele norm)))))))

(addtest (genotype-call-tests) valid-alleles-p/1
  (ensure (every #'valid-alleles-p *possible-allele-pairs*)))

(addtest (genotype-call-tests) make-chromosome-p/1
  (let ((chrs (mapcar (lambda (x)
                        (format nil "~d" x)) (iota 10 1))))
    (dolist (chr chrs)
      (let ((pred (make-chromosome-p *example-bpm* chr #'string=)))
        (ensure (funcall pred chr))
        (ensure (not (funcall pred "NOT A VALID CHROMOSOME")))))))

(addtest (genotype-call-tests) gtc-open/close/1
  (with-open-file (stream (merge-pathnames "data/example_0000.gtc")
                          :element-type 'octet)
    (let ((gtc (gtc-open stream)))
      (ensure gtc)
      (ensure (= 3 (version-of gtc)))
      (ensure (= 9 (length (toc-of gtc))))
      (ensure (gtc-close gtc)))))

(addtest (genotype-call-tests) with-gtc/1
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
     (ensure (= 3 (version-of gtc)))
     (ensure (= 9 (length (toc-of gtc))))))

(addtest (genotype-call-tests) data-field-of/1
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
    (ensure (= 10 (data-field-of gtc :num-snps)))
    (ensure (string= "example_0000" (data-field-of gtc :sample-name)))
    (ensure (string= "example.bpm"
                     (data-field-of gtc :snp-manifest)))
    (ensure (= 1 (length (data-field-of gtc :normalization-xforms))))

    (ensure-condition (invalid-argument-error)
      (data-field-of gtc :invalid-field))
      
    (ensure (= 10 (length (data-field-of gtc :x-intensities))))
    (ensure (= 10 (length (data-field-of gtc :y-intensities))))
    (ensure (= 10 (length (data-field-of gtc :genotypes))))
    (ensure (= 10 (length (data-field-of gtc :basecalls))))
    (ensure (= 10 (length (data-field-of gtc :genotype-scores))))))

(addtest (genotype-call-tests) data-field-of/2
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
    (ensure (every #'xform-p
                   (coerce (data-field-of gtc :normalization-xforms) 'list)))))

(addtest (genotype-call-tests) data-field-of/3
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
    (ensure (every (lambda (x)
                     (and (vectorp x) (integerp (aref x 0))))
                   (mapcar (lambda (field)
                             (data-field-of gtc field))
                           '(:x-intensities :y-intensities))))))

(addtest (genotype-call-tests) data-field-of/4
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
    (ensure (every (lambda (x)
                     (and (vectorp x) (stringp (aref x 0))))
                   (mapcar (lambda (field)
                             (data-field-of gtc field))
                           '(:genotypes :basecalls))))))

(addtest (genotype-call-tests) data-field-of/5
  (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
    (let ((scores (data-field-of gtc :genotype-scores)))
      (ensure (and (vectorp scores) (floatp (aref scores 0)))))))

(addtest (genotype-call-tests) sim-open-close/1
  (with-open-file (stream *example-normalized-sim* :element-type 'octet)
    (let ((sim (sim-open stream)))
      (ensure sim)
      (ensure (= 1 (version-of sim)))
      (ensure (sim-close sim)))))

(addtest (genotype-call-tests) sim-closed-stream/1
  (ensure-condition (invalid-operation-error)
    (make-instance 'sim :stream (with-open-file
                                    (stream *example-normalized-sim*
                                            :element-type 'octet)
                                  stream))))

(addtest (genotype-call-tests) with-sim/1
  (with-sim (sim *example-normalized-sim*)
    (ensure sim)
    (ensure (= 1 (version-of sim)))
    (ensure (= 255 (name-size-of sim)))
    (ensure (= 5 (num-samples-of sim)))
    (ensure (= 10 (num-probes-of sim)))
    (ensure (= 2 (num-channels-of sim)))
    (ensure (eql 'single-float (format-of sim)))))

(addtest (genotype-call-tests) with-sim/2
  (with-sim (sim *example-raw-sim*)
    (ensure sim)
    (ensure (= 1 (version-of sim)))
    (ensure (= 255 (name-size-of sim)))
    (ensure (= 5 (num-samples-of sim)))
    (ensure (= 10 (num-probes-of sim)))
    (ensure (= 2 (num-channels-of sim)))
    (ensure (eql 'uint16 (format-of sim)))))

(addtest (genotype-call-tests) read-intensities/1
  (with-sim (sim *example-normalized-sim*)
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
     do (with-sim (sim *example-normalized-sim*)
          (compare-intensities (subseq *sim-intensities-hexdump* (* 2 i))
                               (read-intensities sim :start i))))
  ;; :end argument only
  (loop
     for i from 0 below 10
     do (with-sim (sim *example-normalized-sim*)
          (compare-intensities (subseq *sim-intensities-hexdump* 0 (* 2 i))
                               (read-intensities sim :end i))))
  ;; both :start and :end arguments
  (loop
     for i from 0 below 8
     for j from 2 below 10
     do (with-sim (sim *example-normalized-sim*)
          (compare-intensities
           (subseq *sim-intensities-hexdump* (* 2 i) (* 2 j))
           (read-intensities sim :start i :end j)))))

(addtest (genotype-call-tests) read-intensities/3
  ;; single intensity channel set
  (with-sim (sim *example-normalized-sim*)
    (ensure (= (num-channels-of sim)
               (length (read-intensities sim :start 0 :end 1)))))
  ;; start < 0
  (with-sim (sim *example-normalized-sim*)
    (ensure-condition (invalid-argument-error)
      (read-intensities sim :start -1)))
  ;; end == start
  (with-sim (sim *example-normalized-sim*)
    (let ((intensities (read-intensities sim :start 0 :end 0)))
      (ensure (vectorp intensities))
      (ensure (zerop (length intensities)))))
   ;; end < start
   (with-sim (sim *example-normalized-sim*)
     (ensure-condition (invalid-argument-error)
       (read-intensities sim :start 0 :end -1))))

(addtest (genotype-call-tests) skip-intensities/1
  (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
    (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                       "data/example_0002.gtc" "data/example_0003.gtc"
                       "data/example_0004.gtc")))
      (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                     :if-does-not-exist :create :format 'single-float)
        (dolist (file gtc-files)
          (with-gtc (gtc (merge-pathnames file))
            (copy-intensities gtc sim *example-bpm* :normalize t)))))
    (dolist (n '(0 1 2 3 4))
      (let ((expected-name (format nil "example_000~d" n)))
        (with-sim (sim tmp-sim)
          (skip-intensities sim n)
          (multiple-value-bind (intensities name)
              (read-intensities sim)
            (ensure (equal expected-name name)
                    :report "expected ~a, but found ~a"
                    :arguments (expected-name name))))))))

(addtest (genotype-call-tests) copy-intensities/gtc/sim/1
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp :tmpdir (merge-pathnames "data") :type "sim")
      (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                         "data/example_0002.gtc" "data/example_0003.gtc"
                         "data/example_0004.gtc")))
          (with-sim (sim tmp :direction :output :if-exists :supersede
                         :if-does-not-exist :create :format 'single-float)
            (dolist (file gtc-files)
              (with-gtc (gtc (merge-pathnames file))
                (ensure (copy-intensities gtc sim *example-bpm* :normalize t))))
            (ensure (= (length gtc-files) (num-samples-of sim))
                    :report "expected ~d, but found ~d"
                    :arguments ((length gtc-files) (num-samples-of sim))))
          (with-sim (sim tmp)
            (ensure (= (length gtc-files) (num-samples-of sim))))))))

(addtest (genotype-call-tests) mismatched-manifest/1
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp :tmpdir (merge-pathnames "data") :type "sim")
      (with-sim (sim tmp :direction :output :if-exists :supersede
                     :if-does-not-exist :create)
        (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
          (ensure-condition (invalid-argument-error)
              (copy-intensities gtc sim *alternative-bpm*)))))))

(addtest (genotype-call-tests) sim-to-illuminus/1
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
      (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                         "data/example_0002.gtc" "data/example_0003.gtc"
                         "data/example_0004.gtc")))
        (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                       :if-does-not-exist :create :format 'single-float)
          (dolist (file gtc-files)
            (with-gtc (gtc (merge-pathnames file))
              (copy-intensities gtc sim *example-bpm* :normalize t))))
        (with-tmp-pathname (tmp-iln :tmpdir (merge-pathnames "data")
                                    :type "iln")
          (let ((iln (sim-to-illuminus tmp-iln *example-bpm* tmp-sim)))
            (ensure iln))
          (ensure-lines-equal
           (merge-pathnames "data/example.iln") tmp-iln))))))

(addtest (genotype-call-tests) sim-to-illuminus/2
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
      (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                     :if-does-not-exist :create :format 'single-float)
        (with-gtc (gtc (merge-pathnames "data/example_0000.gtc"))
          (copy-intensities gtc sim *example-bpm*
                            :test (lambda (index)
                                    (oddp index))
                            :key #'snp-index
                            :normalize t)))
      (with-tmp-pathname (tmp-iln :tmpdir (merge-pathnames "data") :type "iln")
        (ensure (sim-to-illuminus tmp-iln *example-bpm* tmp-sim
                                  :test #'(lambda (index)
                                            (oddp index))
                                  :key #'snp-index))
        (ensure-lines-equal
         (merge-pathnames "data/example_odd.iln") tmp-iln)))))

(addtest (genotype-call-tests) sim-to-illuminus/3
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
      (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                         "data/example_0002.gtc" "data/example_0003.gtc"
                         "data/example_0004.gtc")))
        (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                       :if-does-not-exist :create :format 'single-float)
          (dolist (file gtc-files)
            (with-gtc (gtc (merge-pathnames file))
              (copy-intensities gtc sim *example-bpm* :normalize t))))
        ;; Test restricting to specific chromosomes
        (dolist (chr (chromosomes-of *example-bpm*))
          (multiple-value-bind (start end)
              (chromosome-boundaries *example-bpm* chr)
            (with-tmp-pathname (tmp-iln :tmpdir (merge-pathnames "data")
                                        :type "iln")
              (let ((iln (sim-to-illuminus tmp-iln *example-bpm* tmp-sim
                                           :start start :end end)))
                (ensure iln))
              (let ((observed (with-open-file (s tmp-iln)
                                (loop
                                   for line = (read-line s nil nil)
                                   while line
                                   collect line into lines
                                   finally (return (rest lines)))))
                    (expected (with-open-file (s (merge-pathnames
                                                  "data/example.iln"))
                                (loop
                                   for line = (read-line s nil nil)
                                   while line
                                   collect line into lines
                                   finally (return (subseq (rest lines)
                                                           start end))))))
                (ensure (equalp expected observed)
                        :report "expected ~a, but found ~a"
                        :arguments (expected observed))))))))))

(addtest (genotype-call-tests) sim-to-genosnp/1
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
      (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                         "data/example_0002.gtc" "data/example_0003.gtc"
                         "data/example_0004.gtc")))
        (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                       :if-does-not-exist :create :format 'uint16)
          (dolist (file gtc-files)
            (with-gtc (gtc (merge-pathnames file))
              (copy-intensities gtc sim *example-bpm* :normalize nil))))
        (with-tmp-pathname (tmp-gsn :tmpdir (merge-pathnames "data")
                                    :type "gsn")
          (let ((gsn (sim-to-genosnp tmp-gsn *example-bpm* tmp-sim)))
            (ensure gsn))
          (ensure-lines-equal
           (merge-pathnames "data/example.raw.gsn") tmp-gsn))))))

(addtest (genotype-call-tests) sim-to-genosnp/2
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-sim :tmpdir (merge-pathnames "data") :type "sim")
      (let ((gtc-files '("data/example_0000.gtc" "data/example_0001.gtc"
                         "data/example_0002.gtc" "data/example_0003.gtc"
                         "data/example_0004.gtc")))
        (with-sim (sim tmp-sim :direction :output :if-exists :supersede
                       :if-does-not-exist :create :format 'uint16)
          (dolist (file gtc-files)
            (with-gtc (gtc (merge-pathnames file))
              (copy-intensities gtc sim *example-bpm* :normalize nil))))
        (with-tmp-pathname (tmp-gsn :tmpdir (merge-pathnames "data")
                                    :type "gsn")
          (let ((gsn (sim-to-genosnp tmp-gsn *example-bpm* tmp-sim
                                     :start 1 :end 4)))
            (ensure gsn))
          (ensure-lines-equal
           (merge-pathnames "data/example.raw.skipped.gsn") tmp-gsn))))))

(addtest (genotype-call-tests) gtc-to-bed/1
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-bed :tmpdir (merge-pathnames "data") :type "bed")
      (let ((specs (example-gtc-specs))
            (bim-file (plink-pathname tmp-bed "bim"))
            (fam-file (plink-pathname tmp-bed "fam")))
        (ensure (gtc-to-bed tmp-bed *example-bpm* specs))
        (ensure (binary-file= tmp-bed (merge-pathnames "data/example.bed")))
        (ensure (probe-file bim-file))
        (ensure (probe-file fam-file))
        (ensure-lines-equal bim-file (merge-pathnames "data/example.bim"))
        (ensure-lines-equal fam-file (merge-pathnames "data/example.fam"))
        (delete-file bim-file)
        (delete-file fam-file)))))

(addtest (genotype-call-tests) gtc-to-bed/2
  (handler-bind ((test-condition #'leave-tmp-pathname))
    (with-tmp-pathname (tmp-bed :tmpdir (merge-pathnames "data") :type "bed")
      (let ((specs (example-gtc-specs))
            (bim-file (plink-pathname tmp-bed "bim"))
            (fam-file (plink-pathname tmp-bed "fam")))
        (ensure (gtc-to-bed tmp-bed *example-bpm* specs
                            :test #'(lambda (index)
                                      (oddp index))
                            :key #'snp-index))
        (ensure (binary-file= tmp-bed (merge-pathnames "data/example_odd.bed")))
        (ensure (probe-file bim-file))
        (ensure (probe-file fam-file))
        (ensure-lines-equal bim-file (merge-pathnames "data/example_odd.bim"))
        (ensure-lines-equal fam-file (merge-pathnames "data/example_odd.fam"))
        (delete-file bim-file)
        (delete-file fam-file)))))
