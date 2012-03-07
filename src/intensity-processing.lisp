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

;; Implementation of GTC -> SIM with SNP vector metadata
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

;;; Implementation of GTC -> SIM with SNP manifest metadata
(defmethod copy-intensities ((gtc gtc) (sim sim) (manifest bpm)
                             &key key test)
  (copy-intensities gtc sim (snps-of manifest :key key :test test)))

;;; Implementation of SIM -> Illuminus with SNP vector metadata
(defmethod copy-intensities ((sim sim) (iln iln) (snps vector)
                             &key key test (start 0) end)
  (with-slots (num-samples num-probes num-channels)
      sim
    (let ((snps (if test
                    (remove-if test snps :key key)
                    snps))
          (end (or end num-probes)))
      (check-arguments (and (integerp start) (not (minusp start))) (start)
                       "expected a non-negative integer")
      (check-arguments (and (integerp end) (not (minusp end))) (end)
                       "expected a non-negative integer")
      (check-arguments (<= 0 start end) (start end)
                       "start and end must satisfy 0 <= start <= end")
      (check-arguments (= num-probes (length snps)) (sim snps key test)
                       "~d annotations were selected for ~d probes"
                       (length snps) num-probes)
      (let ((sample-names (make-array num-samples))
            (intensities (make-array num-samples)))
        (loop
           for i from 0 below num-samples
           do (multiple-value-bind (sample-intensities name)
                  (read-intensities sim :start start :end end)
                (setf (svref sample-names i) name
                      (svref intensities i) sample-intensities)))
        ;; Should check here that all the intensity vectors are the
        ;; same length
        (let ((stream (stream-of iln))
              (*print-pretty* nil))
          (write-illuminus-header sample-names stream)
          ;; For each SNP
          (loop
             for i from 0 below end     ; probe index, offset to 0
             for j = (* 2 i)            ; intensity index, in pairs
             for k = start then (1+ k)  ; SNP index, not offset
             do (progn
                  ;; Write this SNP's intensities for all samples
                  (write-illuminus-snp (svref snps k) stream)
                  (loop                
                     for sample-intensities across intensities
                     do (write-illuminus-intensities
                         (aref sample-intensities j)
                         (aref sample-intensities (1+ j))
                         stream)
                     finally (terpri stream))))))))
  iln)

(defun num-intensities (intensities)
  (let ((lengths (map 'list #'length intensities)))
    (cond ((null lengths)
           0)
          ((every (lambda (x)
                    (= (first lengths) x)) lengths)
           (first lengths))
          (t
           (error 'invalid-operation-error
                  :format-control "intensity vectors were different ~
                                   lengths: ~a"
                  :format-arguments (list lengths))))))


;;; Implementation of SIM -> Illuminus with SNP manifest metadata
(defmethod copy-intensities ((sim sim) (iln iln) (manifest bpm)
                             &key key test (start 0) end)
  (copy-intensities sim iln (snps-of manifest :key key :test test)
                    :start start :end end))

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

(defgeneric sim-to-illuminus (illuminus-filespec manifest sim-filespec
                              &key test key start end)
  (:method (illuminus-filespec (manifest bpm) sim-filespec
            &key test key (start 0) end)
    (with-sim (sim sim-filespec)
      (if (streamp illuminus-filespec)
          (copy-intensities sim (make-instance 'iln :stream illuminus-filespec)
                            manifest :key key :test test
                            :start start :end end)
          (with-open-file (stream illuminus-filespec :direction :output
                                  :external-format :ascii
                                  :if-exists :supersede
                                  :if-does-not-exist :create)
        (copy-intensities sim (make-instance 'iln :stream stream) manifest
                          :key key :test test :start start :end end))))))
