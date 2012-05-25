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
                                     &key key test name normalize)
  (declare (ignorable key test name normalize))
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
                             &key key test name normalize)
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
      (check-arguments (or (and normalize (= 4 isize))
                           (and (not normalize) (= 2 isize))) (sim normalize)
                           (txt "SIM must be 4 bytes per channel for"
                                "normalized and 2 bytes per channel for raw"
                                "intensities"))
      (write-sequence (string-to-octets sample-name) stream)
      (dotimes (n (- name-size (length sample-name))) ; pad the name
        (write-byte 0 stream))
      (if normalize
          (write-norm-2channel-intensities snps x-intensities y-intensities
                                           xforms stream)
          (write-raw-2channel-intensities snps x-intensities y-intensities
                                          stream))))
  sim)

(defmethod copy-intensities :after ((gtc gtc) (sim sim) (snps vector)
                                    &key key test name normalize)
  (declare (ignorable key test name normalize))
  (with-slots (num-samples)
      sim
    (incf num-samples)))

;;; Implementation of GTC -> SIM with SNP manifest metadata
(defmethod copy-intensities ((gtc gtc) (sim sim) (manifest bpm)
                             &key key test name normalize)
  (let ((gtc-name (manifest-name-of gtc))
        (man-name (name-of manifest)))
    (when (and gtc-name man-name)
      (check-arguments (string= man-name gtc-name)
                       (gtc manifest)
                       "manifest ~s does not match the expected manifest ~s"
                        gtc-name man-name)))
  (copy-intensities gtc sim (snps-of manifest :key key :test test)
                    :name name :normalize normalize))

;;; Implementation of SIM -> Illuminus with SNP vector metadata
(defmethod copy-intensities ((sim sim) (iln iln) (snps vector)
                             &key (start 0) end)
  (with-slots (num-samples num-probes)
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

;;; Implementation of SIM -> GenoSNP with SNP vector metadata
(defmethod copy-intensities ((sim sim) (gsn gsn) (snps vector)
                             &key (start 0) end)
  (with-slots (num-samples num-probes)
      sim
    (let ((end (or end num-samples)))
      (check-arguments (and (integerp start) (not (minusp start))) (start)
                       "expected a non-negative integer")
      (check-arguments (and (integerp end) (not (minusp end))) (end)
                       "expected a non-negative integer")
      (check-arguments (<= 0 start end) (start end)
                       "start and end must satisfy 0 <= start <= end")
      (check-arguments (= num-probes (length snps)) (sim snps)
                       "SIM holds data for ~d SNPS, but found ~d"
                       num-probes (length snps))
      (let ((stream (stream-of gsn))
            (*print-pretty* nil))
        ;; For each sample
        (loop
           for i from start below end
           do (multiple-value-bind (sample-intensities sample-name)
                  (read-intensities sim)
                (write-genosnp-sample sample-name stream)
                (loop
                   for j from 0 below num-probes
                   for k from 0 by 2    ; intensity index, in pairs
                   do (write-genosnp-intensities
                       (aref sample-intensities k)
                       (aref sample-intensities (1+ k))
                       stream)
                   finally (terpri stream)))))))
  gsn)

;;; Implementation of SIM -> GenoSNP with SNP manifest metadata
(defmethod copy-intensities ((sim sim) (gsn gsn) (manifest bpm)
                             &key key test (start 0) end)
  (copy-intensities sim gsn (snps-of manifest :key key :test test)
                    :start start :end end))


(defgeneric gtc-to-sim (sim-filespec manifest sample-specs
                        &key test key normalize)
  (:documentation "Creates a new SIM file containing the aggregated
intensity data from a list of GTC files. Writes a JSON file containing the
chromsome boundaries (in terms of SNP columns)")
  (:method (sim-filespec (manifest bpm) sample-specs &key test key normalize)
    (with-sim (sim sim-filespec :direction :output :if-exists :supersede
                   :if-does-not-exist :create :format (if normalize
                                                          'single-float
                                                          'uint16))
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
                             (puri:urn-nss uri)
                             "conflict with sample name in GTC file '~a'"
                             gtc-name))
              (copy-intensities gtc sim snps :name (format nil "~a" uri)
                                :normalize normalize)))))
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

(defgeneric sim-to-genosnp (genosnp-filespec manifest sim-filespec
                            &key test key start end)

  (:method (genosnp-filespec (manifest bpm) sim-filespec
            &key test key (start 0) end)
    (with-sim (sim sim-filespec)
      (if (streamp genosnp-filespec)
          (copy-intensities sim (make-instance 'gsn :stream genosnp-filespec)
                            manifest :key key :test test
                            :start start :end end)
          (with-open-file (stream genosnp-filespec :direction :output
                                  :external-format :ascii
                                  :if-exists :supersede
                                  :if-does-not-exist :create)
            (copy-intensities sim (make-instance 'gsn :stream stream) manifest
                              :key key :test test :start start :end end))))))
