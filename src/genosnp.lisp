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

(defun write-genosnp-sample (sample-name stream)
  (write-string sample-name stream)
  (write-char #\Tab stream)
  (write-string sample-name stream)
  (write-char #\Tab stream))

(defun write-genosnp-intensities (a b stream)
  (write-char #\Tab stream)
  (princ a stream)
  (write-char #\Space stream)
  (princ b stream))
