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

(defgeneric copy-intensities (from to metadata &key key test)
  (:documentation ""))

(defmethod copy-intensities :before ((gtc gtc) (sim sim) (manifest bpm)
                                     &key key test)
  (with-slots (stream version name-size num-probes num-samples num-channels)
      sim
    (check-arguments (= 2 num-channels) (sim)
                     "found ~d intensity channels; expected 2 for GTC data"
                     num-channels)
    (let ((num-snps (num-snps-of manifest :key key :test test)))
      (check-arguments (= num-probes num-snps) (manifest sim)
                       "SIM holds data for ~d SNPs, but found ~d"
                       num-probes num-snps))))

(defmethod copy-intensities ((gtc gtc) (sim sim) (manifest bpm)
                             &key key test)
  (with-slots (stream name-size)
      sim
    (let ((sample-name (data-field-of gtc :sample-name))
          (xforms (data-field-of gtc :normalization-xforms))
          (x-intensities (data-field-of gtc :x-intensities))
          (y-intensities (data-field-of gtc :y-intensities))
          (isize (intensity-size-of sim))
          (snps (snps-of manifest :key key :test test)))
      (check-arguments (<= (length sample-name) name-size) (gtc sim)
                       (txt "sample name ~s is too long to fit in SIM file"
                            "with limit of ~d characters")
                       sample-name name-size)
      (write-sequence (string-to-octets sample-name) stream)
      (dotimes (n (- name-size (length sample-name))) ; pad the name
        (write-byte 0 stream))
      (write-2channel-intensities snps x-intensities y-intensities
                                  isize xforms stream))))

(defmethod copy-intensities :after ((gtc gtc) (sim sim) (manifest bpm)
                                    &key key test)
  (declare (ignore key test))
  (with-slots (num-samples)
      sim
    (incf num-samples)))

(defmethod copy-intensities ((sim sim) (manifest bpm) stream
                             &key key test)
  (with-slots (num-samples num-probes num-channels)
      sim
    (let ((snps (snps-of manifest :key key :test test)))
      (check-arguments (= num-probes (length snps)) (sim manifest key test)
                       "~d annotations were selected for ~d probes"
                       (length snps) num-probes)
      (let ((sample-names (make-array num-samples))
            (sample-intensities (make-array (* num-probes num-channels))))
        (loop
           for i from 0 below num-samples
           do (multiple-value-bind (intensities name)
                  (read-intensities sim)
                (setf (svref sample-names i) name
                      (svref sample-intensities i) intensities)))
        (write-illuminus-header sample-names stream)
        (loop
           for i from 0 below num-probes
           for j from 0 by num-channels
           do (progn
                (write-illuminus-snp (svref snps i) stream)
                (loop
                   for k from 0 below num-samples
                   do (let ((intensities (svref sample-intensities k)))
                        (write-illuminus-intensities
                         (aref intensities j) (aref intensities (1+ j))
                         stream)))
                (terpri stream)))))))
