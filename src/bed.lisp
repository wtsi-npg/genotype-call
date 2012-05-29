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

(defvar *bed-magic* (make-array 2 :element-type 'octet
                                :initial-contents '(#b01101100 #b00011011))
  "The Plink BED file magic number.")

(defvar *bed-snp-major* #b00000001
  "The Plink BED file flag byte for SNP-major orientation.")
(defvar *bed-individual-major* #b00000000
  "The Plink BED file flag byte for individual-major orientation.")

(defvar *plink-unknown* "-9")

(defclass bed ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The input or output stream.")
   (orientation :initform :individual-major :initarg :orientation
                :reader orientation-of
                :documentation "The file mode indicating SNP-major or
                individual-major orientation. Defaults to
                individual-major."))
  (:documentation "A Plink binary PED (pedigree) file."))

(defmethod initialize-instance :after ((bed bed) &key)
  (with-slots (stream orientation)
      bed
    (cond ((input-stream-p stream)
           (read-magic stream (make-array 2 :element-type 'octet
                                          :initial-element 0) *bed-magic*)
           (setf orientation (read-bed-orientation stream)))
          ((output-stream-p stream)
           (write-magic *bed-magic* stream)
           (write-bed-orientation orientation stream)))))

(defmethod print-object ((bed bed) stream)
  (print-unreadable-object (bed stream :type t :identity nil)
    (with-accessors ((s stream-of) (orientation orientation-of))
        bed
      (format stream "~@[~a ~][~a orientation]"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s)) orientation))))

(defmacro with-bed ((var filespec &rest args) &body body)
  `(let ((,var (bed-open ,filespec ,@args)))
     (unwind-protect
          (progn
            ,@body)
       (when ,var
         (bed-close ,var)))))

(defun bed-open (filespec &rest args)
  (multiple-value-bind (args initargs)
      (remove-key-values '(:orientation) args)
    (let ((stream (apply #'open filespec :element-type 'octet args)))
      (if (open-stream-p stream)
          (apply #'make-instance 'bed :stream stream (flatten initargs))
          (error 'invalid-operation-error
                 :format-control "BED stream ~a closed unexpectedly"
                 :format-arguments (list stream))))))

(defun bed-close (bed)
  (close (stream-of bed)))


(defun read-bed-orientation (stream)
  (let ((orientation (read-byte stream nil nil)))
    (check-record (and (numberp orientation)
                       (or (= *bed-snp-major* orientation)
                           (= *bed-individual-major* orientation)))
                  "invalid BED orientation; expected ~a or ~a"
                  *bed-snp-major* *bed-individual-major*)
    (if (= *bed-snp-major* orientation)
        :snp-major
        :individual-major)))

(defun write-bed-orientation (orientation stream)
  (write-byte (ecase orientation
                (:snp-major *bed-snp-major*)
                (:individual-major *bed-individual-major*)) stream))

(defun read-bed-genotypes (stream)
  (flet ((decode-genotype (x)
           (ecase x
             (#b00 "AA")
             (#b11 "BB")
             (#b10 "AB")                ; read LS bit first i.e. 01
             (#b01 nil))))              ; read LS bit first i.e. 10
    (let ((byte (read-byte stream t)))
      (loop
         for dibit in '(0 2 4 6)
         collect (decode-genotype (ldb (byte 2 dibit) byte))))))

(defun write-bed-genotypes (genotypes stream)
  "Writes GENOTYPES to an octet STREAM, packing 4 genotypes per
byte. The terminal byte may be partially filled, with any unused space
being occupied by 0 (see the Plink BED format specification).

Arguments:

- genotypes (simple-vector string): The genotypes for an individual or SNP.
- stream (output stream): The octet stream.

Returns:
- The number of bytes written."
  (flet ((encode-genotype (x)
           (cond ((null x)
                  #b01)
                 ((string= "AA" x)
                  #b00)
                 ((string= "BB" x)
                  #b11)
                 ((or (string= "AB" x) (string= "BA" x))
                  #b10)
                 (t
                  (check-record nil x "invalid genotype")))))
    (let* ((num-genotypes (length genotypes))
           (bytes (make-array (encoded-bed-size num-genotypes)
                              :element-type 'octet :initial-element 0)))
      (loop
         for i from 0 below num-genotypes by 4
         for j = 0 then (1+ j)
         do (loop
               with byte = 0
               for dibit in '(0 2 4 6)
               for k from i below num-genotypes
               do (setf (ldb (byte 2 dibit) byte)
                        (encode-genotype (svref genotypes k)))
               finally (setf (aref bytes j) byte))
         finally (return (write-sequence bytes stream))))))

(defun encoded-bed-size (n)
  "Returns the number of whole bytes required to encode N genotypes."
  (ceiling n 4))

(defun encode-bim-chromosome (chr)
  "Returns a string that is the Plink BIM encoding for CHR (affects
heterosomes and mitochondrial)."
  (cond ((string= "X" chr)
         "23")
        ((string= "Y" chr)
         "24")
        ((string= "XY" chr)
         "25")
        ((string= "MT" chr)
         "26")
        (t
         chr)))

(defun write-bim-snp (snp stream)
  "Writes a Plink BIM (BInary Map) record of SNP to STREAM."
  (write-string (encode-bim-chromosome (snp-chromosome snp)) stream)
  (write-char #\Tab stream)
  (write-string (snp-name snp) stream)
  (write-char #\Tab stream)
  (princ 0 stream)                      ; genetic position
  (write-char #\Tab stream)
  (princ (snp-position snp) stream)     ; physical position
  (write-char #\Tab stream)
  (write-char (snp-allele-a snp) stream)
  (write-char #\Tab stream)
  (write-char (snp-allele-b snp) stream)
  (terpri stream))

;; Plink FAM format is the first 6 columns of Plink PED format. Note
;; that it is not possible to start any family identifier with a '#'
;; character because this is used to indicate a comment. (Plink format
;; has no concept of escaping characters).

(defun write-fam-individual (name stream &optional gender)
  "Writes a Plink FAM record of an individual with NAME to STREAM."
  (let ((gender (or gender *plink-unknown*)))
    (write-string name stream)            ; family identifier
    (write-char #\Tab stream)
    (write-string name stream)            ; individual identifier
    (write-char #\Tab stream)
    (write-string *plink-unknown* stream) ; paternal identifier
    (write-char #\Tab stream)
    (write-string *plink-unknown* stream) ; maternal identifier
    (write-char #\Tab stream)
    (princ gender stream)     ; gender male=1, female=2, unknown=other
    (write-char #\Tab stream)
    (write-string *plink-unknown* stream) ; phenotype
    (terpri stream)))
