;;;
;;; Copyright (c) 2012 Genome Research Ltd. All rights reserved.
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

(in-package :uk.ac.sanger.genotype-call)

(defun generate-sim (filespec manifest &rest sample-names)
  "Writes a fake SIM file. The expected intensity values are
incremented by (n = 1.0) per SNP within a chromosome and by (n = the
number of SNPs in the manifest) between chromosomes. This makes it
possible to predict the correct intensity for each SNP in test data,
given a manifest. The normalization tranform is chosen to be a no-op."
  (let* ((snps (snps-of manifest))
         (name-size (reduce #'max (mapcar #'length sample-names)))
         (xform (make-xform 1 0.f0 0.f0 1.f0 1.f0)) ; no-op xform
         (xforms (make-array 1 :initial-element xform))
         (n (num-snps-of manifest)))
    (with-sim (sim filespec :direction :output :if-exists :supersede
                   :if-does-not-exist :create
                   :num-samples (length sample-names)
                   :num-probes (length snps)
                   :name-size name-size)
      (let ((stream (stream-of sim))
            (isize (intensity-size-of sim))
            (x-intensities (make-array n :element-type 'uint16
                                       :initial-element 0))
            (y-intensities (make-array n :element-type 'uint16
                                       :initial-element 0)))
        (loop
           with snps-per-chr = (make-hash-table :test #'equal)
           for snp across snps
           for i from 0
           do (let* ((chr (snp-chromosome snp))
                     (ci (position chr (chromosomes-of manifest)
                                   :test #'string=))
                     (si (1+ (gethash chr snps-per-chr 0)))
                     (ii (+ (* n ci) si)))
                (setf (gethash chr snps-per-chr) si
                      (aref x-intensities i) ii
                      (aref y-intensities i) ii)))
        (dolist (sample-name sample-names)
          (write-sequence (string-to-octets sample-name) stream)
          (dotimes (len (- name-size (length sample-name))) ; pad the name
            (write-byte 0 stream))
          (write-2channel-intensities
           snps x-intensities y-intensities isize xforms stream))))))
