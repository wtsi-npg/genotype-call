;;;
;;; Copyright (c) 2011 Genome Research Ltd. All rights reserved.
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

(defparameter *example-snp-counts*
  '(("16" 4) ("2" 4) ("9" 2)))

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
                (ensure (write-intensities gtc manifest sim)))
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

(addtest (genotype-call-tests) with-sim/1
  (with-sim (sim (merge-pathnames "data/example.sim"))
    (ensure sim)
    (ensure (= 1 (version-of sim)))
    (ensure (= 255 (name-size-of sim)))
    (ensure (= 5 (num-samples-of sim)))
    (ensure (= 10 (num-probes-of sim)))
    (ensure (= 2 (num-channels-of sim)))
    (ensure (eql 'single-float (format-of sim)))))
