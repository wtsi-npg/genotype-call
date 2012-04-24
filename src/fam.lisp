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

(defvar *plink-unknown* "-9")

;; Plink FAM format is the first 6 columns of Plink PED format. Note
;; that it is not possible to start any family identifier with a '#'
;; character because this is used to indicate a comment. (Plink format
;; has no concept of escaping characters).

(defun write-fam-individual (name stream)
  "Writes a Plink FAM record of an individual with NAME to STREAM."
  (write-string name stream)            ; family identifier
  (write-char #\Tab stream)
  (write-string name stream)            ; individual identifier
  (write-char #\Tab stream)
  (write-string *plink-unknown* stream) ; paternal identifier
  (write-char #\Tab stream)
  (write-string *plink-unknown* stream) ; maternal identifier
  (write-char #\Tab stream)
  (write-string *plink-unknown* stream) ; gender male=1, female=2, unknown=other
  (write-char #\Tab stream)
  (write-string *plink-unknown* stream) ; phenotype
  (terpri stream))                      
