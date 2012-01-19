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

(defgeneric copy-intensities (from to metadata &key test key)
  (:documentation "Copies microarray intensity data FROM sources such
as GTC or SIM files TO a destination such as a SIM or genotype caller
input file. The METADATA argument is used to provide additional
information required during the operation.

Arguments:

- from (object): A producer of intensity data, such as a GTC object
  SIM object or a Lisp stream.
- to (object): A consumer of intensity data, such as a SIM object or
  Lisp stream.
- metadata (object): Metadata used to identify SNPs and samples, for
  example a BeadPool Manifest (BPM object).

Key:

- test (predicate): A test predicate used against the metadata to
  select intensities for inclusion in the output.
- key (function): A function used to transform metadata elements
  before applying TEST.

Returns:

- to (object)."))

(defmethod copy-intensities :before ((gtc gtc) (sim sim) (snps vector)
                                     &key key test)
  (declare (ignorable key test))
  (with-slots (stream version name-size num-probes num-samples num-channels)
      sim
    (check-arguments (= 2 num-channels) (sim)
                     "found ~d intensity channels; expected 2 for GTC data"
                     num-channels)
    (let ((num-snps (length snps)))
      (if (zerop num-samples)
          (setf num-probes num-snps)
          (check-arguments (= num-probes num-snps) (sim snps)
                           "SIM holds data for ~d SNPs, but found ~d"
                           num-probes num-snps)))))

(defmethod copy-intensities ((gtc gtc) (sim sim) (snps vector)
                             &key key test)
  (with-slots (stream name-size)
      sim
    (let ((sample-name (data-field-of gtc :sample-name))
          (xforms (data-field-of gtc :normalization-xforms))
          (x-intensities (data-field-of gtc :x-intensities))
          (y-intensities (data-field-of gtc :y-intensities))
          (isize (intensity-size-of sim))
          (snps (if test
                    (remove-if test snps :key key)
                    snps)))
      (check-arguments (<= (length sample-name) name-size) (gtc sim)
                       (txt "sample name ~s is too long to fit in SIM file"
                            "with limit of ~d characters")
                       sample-name name-size)
      (write-sequence (string-to-octets sample-name) stream)
      (dotimes (n (- name-size (length sample-name))) ; pad the name
        (write-byte 0 stream))
      (write-2channel-intensities snps x-intensities y-intensities
                                  isize xforms stream)))
  sim)

(defmethod copy-intensities :after ((gtc gtc) (sim sim) (snps vector)
                                    &key key test)
  (declare (ignorable key test))
  (with-slots (num-samples)
      sim
    (incf num-samples)))

(defmethod copy-intensities ((gtc gtc) (sim sim) (manifest bpm)
                             &key key test)
  (copy-intensities gtc sim (snps-of manifest :key key :test test)))

(defmethod copy-intensities ((sim sim) stream (manifest bpm)
                             &key key test (start 0) end)
  (with-slots (num-samples num-probes num-channels)
      sim
    (let ((snps (snps-of manifest :key key :test test)))
      (check-arguments (= num-probes (length snps)) (sim manifest key test)
                       "~d annotations were selected for ~d probes"
                       (length snps) num-probes)
      (let ((sample-names (make-array num-samples))
            (intensities (make-array num-samples)))
        (loop
           for i from 0 below num-samples
           do (multiple-value-bind (sample-intensities name)
                  (read-intensities sim :start start :end end)
                (princ ".")
                (setf (svref sample-names i) name
                      (svref intensities i) sample-intensities)))
        (let ((*print-pretty* nil))
          (write-illuminus-header sample-names stream)
          (loop
             for i from start below end by num-channels
             for j = (1+ i)
             do (progn
                  (write-illuminus-snp (svref snps i) stream)
                  (loop
                     for sample-intensities across intensities
                     do (write-illuminus-intensities
                         (aref sample-intensities i) (aref sample-intensities j)
                         stream)
                     finally (terpri stream))))))))
  stream)

(defgeneric gtc-to-sim (sim-filespec manifest gtc-filespecs &key test key)
  (:documentation "Creates a new SIM file containing the aggregated
intensity data from a list of GTC files.")
  (:method (sim-filespec (manifest bpm) gtc-filespecs &key test key)
    (with-sim (sim sim-filespec :direction :output :if-exists :supersede
                   :if-does-not-exist :create)
      (let ((snps (snps-of manifest :test test :key key))
            (prev-manifest-name))
        (dolist (gtc-filespec gtc-filespecs sim)
          (with-gtc (gtc gtc-filespec)
            (let ((manifest-name (data-field-of gtc :snp-manifest)))
              (if (null prev-manifest-name)
                  (setf prev-manifest-name manifest-name)
                  (check-arguments (string= prev-manifest-name manifest-name)
                                   (gtc-filespecs)
                                   "manifest ~s does not match previous ~
                                  manifest ~s" manifest-name
                                  prev-manifest-name)))
            (copy-intensities gtc sim snps))))
      sim)))

;; (defgeneric sim-to-illuminus (illuminus-filespec manifest sim-filespecs
;;                                                  &key test key)
;;   )
