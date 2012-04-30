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

(defgeneric copy-intensities (from to metadata &key key test &allow-other-keys)
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

- key (function): A function used to transform metadata elements
  before applying TEST.
- test (predicate): A test predicate used against the metadata to
  select intensities for inclusion in the output.

Returns:

- to (object)."))

;; Implementation of GTC -> SIM with SNP vector metadata
(defmethod copy-intensities :before ((gtc gtc) (sim sim) (snps vector)
                                     &key key test name)
  (declare (ignorable key test name))
  (with-slots (stream version name-size num-samples num-probes num-channels)
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
                             &key key test name)
  (with-slots (stream name-size)
      sim
    (let ((sample-name (or name (data-field-of gtc :sample-name)))
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
                                    &key key test name)
  (declare (ignorable key test name))
  (with-slots (num-samples)
      sim
    (incf num-samples)))

;;; Implementation of GTC -> SIM with SNP manifest metadata
(defmethod copy-intensities ((gtc gtc) (sim sim) (manifest bpm)
                             &key key test name)
  (let ((gtc-name (manifest-name-of gtc))
        (man-name (name-of manifest)))
    (when (and gtc-name man-name)
      (check-arguments (string= man-name gtc-name)
                       (gtc manifest)
                       "manifest ~s does not match the expected manifest ~s"
                        gtc-name man-name)))
  (copy-intensities gtc sim (snps-of manifest :key key :test test) :name name))

;;; Implementation of SIM -> Illuminus with SNP vector metadata
(defmethod copy-intensities ((sim sim) (iln iln) (snps vector)
                             &key (start 0) end)
  (with-slots (num-samples num-probes num-channels)
      sim
    (let ((end (or end num-probes)))
      (check-arguments (and (integerp start) (not (minusp start))) (start)
                       "expected a non-negative integer")
      (check-arguments (and (integerp end) (not (minusp end))) (end)
                       "expected a non-negative integer")
      (check-arguments (<= 0 start end) (start end)
                       "start and end must satisfy 0 <= start <= end")
      (check-arguments (= num-probes (length snps)) (sim snps)
                       "SIM holds data for ~d SNPS, but found ~d"
                       num-probes (length snps))
      (let ((sample-names (make-array num-samples))
            (intensities (make-array num-samples)))
        (loop
           for i from 0 below num-samples
           do (multiple-value-bind (sample-intensities name)
                  (read-intensities sim :start start :end end)
                (setf (svref sample-names i) name
                      (svref intensities i) sample-intensities)))
        (let ((stream (stream-of iln))
              (*print-pretty* nil))
          (write-illuminus-header sample-names stream)
          ;; For each SNP
          (loop
             for i from start below end ; SNP index, offset to start
             for j from 0 by 2          ; intensity index, in pairs, not offset
             do (progn
                  ;; Write this SNP's intensities for all samples
                  (write-illuminus-snp (svref snps i) stream)
                  (loop                
                     for sample-intensities across intensities
                     do (write-illuminus-intensities
                         (aref sample-intensities j)
                         (aref sample-intensities (1+ j))
                         stream)
                     finally (terpri stream))))))))
  iln)


;;; Implementation of SIM -> Illuminus with SNP manifest metadata
(defmethod copy-intensities ((sim sim) (iln iln) (manifest bpm)
                             &key key test (start 0) end)
  (copy-intensities sim iln (snps-of manifest :key key :test test)
                    :start start :end end))

(defgeneric gtc-to-sim (sim-filespec manifest sample-specs &key test key)
  (:documentation "Creates a new SIM file containing the aggregated
intensity data from a list of GTC files.")
  (:method (sim-filespec (manifest bpm) sample-specs &key test key)
    (with-sim (sim sim-filespec :direction :output :if-exists :supersede
                   :if-does-not-exist :create)
      (let ((snps (snps-of manifest :test test :key key))
            (prev-manifest-name (name-of manifest)))
        (dolist (spec sample-specs sim)
          (let* ((gtc-filespec (assocdr :result spec))
                 (uri (assocdr :uri spec)))
            (with-gtc (gtc gtc-filespec)
              (let ((manifest-name (data-field-of gtc :snp-manifest))
                    (gtc-name (data-field-of gtc :sample-name)))
                (if (null prev-manifest-name)
                    (setf prev-manifest-name manifest-name)
                    (check-arguments (string= prev-manifest-name manifest-name)
                                     (spec)
                                     "manifest ~s does not match the expected ~
                                      manifest ~s" manifest-name
                                      prev-manifest-name))
                (check-field (or (string= "" gtc-name)
                                 (string= (puri:urn-nss uri) gtc-name))
                             uri
                             "conflict with sample name in GTC file"))
              (copy-intensities gtc sim snps :name (format nil "~a" uri))))))
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
