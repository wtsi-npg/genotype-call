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
   (version :initform 1 :reader version-of
            :documentation "The SIM format version.")
   (name-size :initform 255 :initarg :name-size :reader name-size-of)
   (num-probes :initform 0 :initarg :num-probes :reader num-probes-of)
   (num-samples :initform 0 :initarg :num-samples :accessor num-samples-of)
   (num-channels :initform 2 :initarg :num-channels :reader num-channels-of)
   (format :initform 'single-float :initarg :format :reader format-of)))

(defmethod initialize-instance :after ((sim sim) &key)
  (with-slots (stream version name-size num-probes num-samples
                      num-channels format)
      sim
    (cond ((not (open-stream-p stream))
           (error 'invalid-operation-error
                  :format-control "SIM stream ~a  closed unexpectedly"
                  :format-arguments (list stream)))
          ((input-stream-p stream)
           (setf (values version name-size num-probes num-samples
                         num-channels format)
                 (read-sim-header stream (make-array 4 :element-type 'octet
                                                     :initial-element 0))))
          ((output-stream-p stream)
           (write-sim-header version name-size num-probes num-samples
                             num-channels format stream
                             (make-array 4 :element-type 'octet))))))

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
  (multiple-value-bind (oargs initargs)
      (remove-key-values '(:version :name-size :num-probes
                           :num-samples :num-channels :format) args)
    (let ((stream (apply #'open filespec :element-type 'octet oargs)))
      (apply #'make-instance 'sim :stream stream (flatten initargs)))))

(defun sim-close (sim)
  (with-slots (stream version name-size num-probes num-samples
                      num-channels format)
      sim
    (when (output-stream-p stream)
      (unless (file-position stream 0)
        (error 'io-error "failed to write header"))
      (write-sim-header version name-size num-probes num-samples
                        num-channels format stream
                        (make-array 4 :element-type 'octet :initial-element 0)))
    (close stream)))

(defgeneric header-size-of (file)
  (:method ((sim sim))
    ;; magic + version + name-size + num-samples + num-probes +
    ;; num-channels + format
    (+ 3 1 2 4 4 1 1)))

(defgeneric intensity-size-of (file)
  (:method ((sim sim))
    (ecase (format-of sim)
      (single-float 4)
      (uint16-scaled 2))))

(defgeneric read-intensities (file &key start end)
  (:method ((sim sim) &key (start 0) end)
    (with-slots (stream name-size num-samples num-probes num-channels)
        sim
      (let ((end (or end num-probes))
            (intensity-size (intensity-size-of sim)))
        (check-arguments (<= start end num-probes) (start end num-probes)
                         "must satisfy start <= end <= number of probes")
        (let ((name-buffer (make-array name-size :element-type 'octet
                                       :initial-element 0))
              (buffer (make-array (* (- end start) num-channels intensity-size)
                                  :element-type 'octet :initial-element 0))
              (intensities (make-array (* (- end start) num-channels)
                                       :element-type 'single-float
                                       :initial-element 0.f0)))
          (unless (= name-size (read-sequence name-buffer stream))
            (error 'malformed-record-error :record name-buffer
                   :format-control "failed to read sample name"))
          ;; Maybe seek to start of intensities
          (when (plusp start)
            (unless (file-position stream (+ (* intensity-size start)
                                             (file-position stream)))
              (error 'malformed-file-error :file sim
                     "failed to seek to start of intensity data")))
          (let ((n (read-sequence buffer stream)))
            (unless (= (length buffer) n)
              (error 'malformed-file-error :file sim
                     :format-control "expected ~d bytes, but read ~d"
                     :format-arguments (list (length buffer) n)))
            (loop
               for i from 0 below n by intensity-size
               for j = 0 then (1+ j)
               do (setf (aref intensities j) (decode-float32le buffer i))))
          ;; Maybe seek to end of intensities
          (when (/= num-probes end)
            (unless (file-position stream (+ (* intensity-size end)
                                             (file-position stream)))
              (error 'malformed-file-error :file sim
                     "failed to seek to end of intensity data")))
          (values intensities (string-right-trim
                               '(#\Nul) (octets-to-string name-buffer))))))))

(defun read-sim-format (stream buffer)
  (ecase (read-uint8 stream buffer)
    (0 'single-float)
    (1 'uint16-scaled)))

(defun write-sim-format (format stream buffer)
  (write-uint8 (ecase format
                 (single-float 0)
                 (uint16-scaled 1)) stream buffer))

(defun read-sim-header (stream buffer)
  (read-magic stream buffer *sim-magic*)
  (values (read-version stream buffer)
          (read-uint16 stream buffer)   ; name-size
          (read-uint32 stream buffer)   ; num-samples
          (read-uint32 stream buffer)   ; num-probes
          (read-uint8 stream buffer)    ; num-channels
          (read-sim-format stream buffer)))

(defun write-sim-header (version name-size num-probes num-samples
                         num-channels format stream buffer)
  (write-magic *sim-magic* stream)
  (write-version version stream buffer)
  (write-uint16 name-size stream buffer)
  (write-uint32 num-probes stream buffer)
  (write-uint32 num-samples stream buffer)
  (write-uint8 num-channels stream buffer)
  (write-sim-format format stream buffer))

(defun write-2channel-intensities (snps x-intensities y-intensities
                                   intensity-size xforms stream)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type simple-vector snps)
           (type (simple-array uint16 (*)) x-intensities y-intensities)
           (type (integer 2 4) intensity-size))
  (let ((total-size (* 2 intensity-size (length snps)))
        (pair-size (* 2 intensity-size)))
    (loop
       with bindata = (make-array total-size :element-type 'octet)
       for snp across snps
       for j of-type fixnum from 0 by pair-size
       for k of-type fixnum = (+ j intensity-size)
       do (let ((m (1- (snp-index snp))))
            (multiple-value-bind (xnorm ynorm)
                (normalize (aref x-intensities m) (aref y-intensities m)
                           (svref xforms (snp-norm-rank snp)))              
              (encode-float32le xnorm bindata j)
              (encode-float32le ynorm bindata k)))
       finally (return (write-sequence bindata stream)))))
