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

(defmacro with-restored-position (stream &body body)
  "Execute BODY, ensuring that the file position of STREAM, if
changed, is restored on leaving."
  (with-gensyms (pos)
    `(let ((,pos (file-position stream)))
       (unwind-protect
            (progn
              ,@body)
         (when (and ,stream (open-stream-p ,stream))
           (unless (= (file-position stream) ,pos)
             (file-position stream ,pos)))))))

(defun read-magic (stream buffer expected)
  (let* ((len (length expected))
         (magic (subseq (read-record stream buffer len) 0  len)))
    (or (equalp expected magic)
        (error 'malformed-file-error :file stream
               :format-control "invalid magic number: expected ~a, found ~a"
               :format-arguments (list expected magic)))))

(defun write-magic (magic stream)
  (write-sequence magic stream))

(defun read-version (stream buffer)
  "Reads the file version from STREAM."
  (decode-uint8le (read-record stream buffer 1)))

(defun write-version (version stream buffer)
  (write-sequence (encode-int8le version buffer) stream :end 1))

(defun read-float (stream buffer)
  (decode-float32le (read-record stream buffer 4)))

(defun read-uint8 (stream buffer)
  (decode-uint8le (read-record stream buffer 1)))

(defun write-uint8 (value stream buffer)
  (write-sequence (encode-int8le value buffer) stream :end 1))

(defun read-uint16 (stream buffer)
  (decode-uint16le (read-record stream buffer 2)))

(defun write-uint16 (value stream buffer)
  (write-sequence (encode-int16le value buffer) stream :end 2))

(defun read-uint32 (stream buffer)
  (decode-uint32le (read-record stream buffer 4)))

(defun write-uint32 (value stream buffer)
  (write-sequence (encode-int32le value buffer) stream :end 4))

(defun read-string (stream buffer)
  (let ((buffer (read-record stream buffer 1)))
    (let* ((len (decode-uint8le buffer))
           (str (make-array len :element-type 'octet :initial-element 0)))
      (octets-to-string (read-record stream str len)))))

(defun read-record (stream buffer record-size)
  (let ((num-bytes (read-sequence buffer stream :end record-size)))
    (unless (= num-bytes record-size)
      (error 'malformed-record-error
             :record (subseq buffer 0 num-bytes)
             :format-control "only read ~d of ~d expected bytes"
             :format-arguments (list num-bytes record-size)))
    buffer))
