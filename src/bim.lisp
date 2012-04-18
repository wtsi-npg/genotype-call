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

(defun write-bim-snp (snp stream)
  "Writes a Plink BIM (BInary Map) record of SNP to STREAM."
  (write-string (snp-chromosome snp) stream)
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
