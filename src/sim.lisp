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

(defvar *sim-magic* (map-into (make-array 3 :element-type 'octet)
                              #'char-code "sim")
  "Simple Intensity matrix (SIM) file magic string.")

(defclass sim ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The stream.")
   (headerp :initform nil)
   (version :initform 1 :reader version-of
            :documentation "The SIM format version.")
   (name-size :initform 255 :initarg :name-size :reader name-size-of)
   (num-samples :initform 0 :initarg :num-samples :accessor num-samples-of)
   (num-probes :initform 0 :initarg :num-probes :accessor num-probes-of)
   (num-channels :initform 2 :initarg :num-channels :reader num-channels-of)
   (format :initform 'single-float :initarg :format :accessor format-of)))

(defmethod initialize-instance :after ((sim sim) &key)
  (with-slots (stream version num-probes num-samples name-size
                      num-channels format)
      sim
    (cond ((not (open-stream-p stream))
           (error 'invalid-operation-error
                  :format-control "SIM stream ~a  closed unexpectedly"
                  :format-arguments (list stream)))
          ((input-stream-p stream)
           (let ((buffer (make-array 4 :element-type 'octet
                                     :initial-element 0)))
             (read-magic stream buffer *sim-magic*)
             (setf version (read-version stream buffer)
                   name-size (read-uint16 stream buffer)
                   num-probes (read-uint32 stream buffer)
                   num-samples (read-uint32 stream buffer)
                   num-channels (read-uint8 stream buffer)
                   format (read-sim-format stream buffer))))
          ((output-stream-p stream)
           t))))

(defmethod print-object ((sim sim) stream)
  (print-unreadable-object (sim stream :type t :identity nil)
    (with-slots ((s stream) version name-size num-samples num-probes
                 num-channels format)
        sim
      (format stream "~@[~a ~][SIM version ~d] ~a, ~d probes/~d channels ~
                      for ~d samples"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s))
              version format num-probes num-channels num-samples))))

(defmacro with-sim ((var filespec &rest args) &body body)
  `(let ((,var (sim-open ,filespec ,@args)))
     (unwind-protect
          (progn
            ,@body)
       (when ,var
         (sim-close ,var)))))

(defun sim-open (filespec &rest args)
  (let ((stream (apply #'open filespec :element-type 'octet args)))
    (make-instance 'sim :stream stream)))

(defun sim-close (sim)
  (with-slots (stream num-samples)
      sim
    (when (output-stream-p stream)
      (unless (file-position stream 10)
        (error 'io-error "failed to write sample count"))
      (write-sequence (encode-int32le
                       num-samples (make-array 4 :element-type 'octet)) stream))
    (close stream)))

(defmethod write-header ((sim sim))
  (with-slots (stream version num-probes num-samples name-size
                      num-channels format headerp)
      sim
    (let ((buffer (make-array 4 :element-type 'octet :initial-element 0)))
      (write-magic *sim-magic* stream)
      (write-version version stream buffer)
      (write-uint16 name-size stream buffer)
      (write-uint32 num-probes stream buffer)
      (write-uint32 num-samples stream buffer)
      (write-uint8 num-channels stream buffer)
      (write-sim-format format stream buffer)
      (setf headerp t))))

(defmethod write-intensities :before ((gtc gtc) (manifest bpm) (sim sim)
                                      &optional chromosome)
  (with-slots (num-probes num-channels headerp)
      sim
    (check-arguments (= 2 num-channels) (sim)
                     "found ~d intensity channels; expected 2 for GTC data"
                     num-channels)
    (let ((num-snps (num-snps-of manifest chromosome)))
      (cond (headerp
             (check-arguments (= num-probes num-snps) (manifest sim)
                              "SIM can hold data for ~d SNPs, but found ~d"
                              num-probes num-snps))
            (t
             (setf num-probes num-snps
                   headerp t)
             (write-header sim))))))

(defmethod write-intensities ((gtc gtc) (manifest bpm) (sim sim)
                              &optional chromosome)
  (let ((stream (stream-of sim))
        (name-len (name-size-of sim))
        (sample-name (data-field-of gtc :sample-name))
        (xforms (data-field-of gtc :normalization-xforms))
        (x-intensities (data-field-of gtc :x-intensities))
        (y-intensities (data-field-of gtc :y-intensities))
        (snps (snps-of manifest chromosome)))
    (check-arguments (<= (length sample-name) name-len) (gtc sim)
                     (txt "sample name ~s is too long to fit in SIM file"
                          "with limit of ~d characters")
                     sample-name name-len)
    (write-sequence (string-to-octets sample-name) stream)
    (dotimes (n (- (name-size-of sim) (length sample-name))) ; pad name
      (write-byte 0 stream))
    (%write-intensities snps x-intensities y-intensities xforms stream)))

(defmethod write-intensities :after ((gtc gtc) (manifest bpm) (sim sim)
                                     &optional chromosome)
  (declare (ignore chromosome))
  (with-slots (num-samples)
      sim
    (incf num-samples)))

(defun read-sim-format (stream buffer)
  (ecase (read-uint8 stream buffer)
    (0 'single-float)
    (1 'uint16-scaled)))

(defun write-sim-format (format stream buffer)
  (write-uint8 (ecase format
                 (single-float 0)
                 (uint16-scaled 1)) stream buffer))

(defun %write-intensities (snps x-intensities y-intensities xforms stream)
  (let* ((len (length snps))
         (bindata (make-array (* 2 4 len) :element-type 'octet)))
    (loop
       for snp across snps
       for j from 0 by 8
       for k = (+ j 4)
       do (let ((m (1- (snp-index snp))))
            (multiple-value-bind (xnorm ynorm)
                (normalize (aref x-intensities m) (aref y-intensities m)
                           (svref xforms (snp-norm-rank snp)))              
              (encode-float32le xnorm bindata j)
              (encode-float32le ynorm bindata k)))
       finally (return (write-sequence bindata stream)))))

