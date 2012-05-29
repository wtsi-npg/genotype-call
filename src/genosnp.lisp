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

;;; GenoSNP intensity (sample) input format
;;; ID1<tab>ID2<tab>1A<space>1B<tab>2A<space>2B<tab> ... nA<space>nB

(defclass gsn ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The GenoSNP intensities stream."))
  (:documentation "GenoSNP format intensities."))

(defmethod print-object ((gsn gsn) stream)
  (print-unreadable-object (gsn stream :type t :identity nil)
    (with-slots ((s stream))
        gsn
      (format stream "~@[~a ~]GenoSNP intensities"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s))))))

(defgeneric save-genosnp-snps (filespec manifest)
  (:documentation "Writes SNP annotation to FILESPEC in the format
expected by GenoSNP. This is a separate method rather than being part
of intensity format-shifting because it is typically done once for all
samples; therefore it is not required most of the time. All SNPs in
the manifest are written.")
  (:method (filespec (manifest bpm))
     (with-open-file (out filespec :direction :output :if-exists :supersede
                          :if-does-not-exist :create)
        (loop
           for snp across (snps-of manifest)
           do (write-genosnp-snp snp out)
             (terpri out)))))

(defun write-genosnp-sample (sample-name stream)
  (write-string sample-name stream)
  (write-char #\Tab stream)
  (write-string sample-name stream)
  (write-char #\Tab stream))

(defun write-genosnp-snp (snp stream)
  (write-string (snp-name snp) stream)
  (write-char #\Tab stream)
  ;; FIXME: The modulo 100 transform here is to emulate the existing
  ;; software (g2i). However, I've no idea why this is applied. In
  ;; fact, GenoSNP docs say that this value should be the BeadPool
  ;; number instead.
  (princ (1+ (mod (snp-norm-id snp) 100)) stream)
  (write-char #\Tab stream)
  (write-char (snp-allele-a snp) stream)
  (write-char #\Space stream)
  (write-char (snp-allele-b snp) stream))

(defun write-genosnp-intensities (a b stream)
  (write-char #\Tab stream)
  (princ a stream)
  (write-char #\Space stream)
  (princ b stream))
