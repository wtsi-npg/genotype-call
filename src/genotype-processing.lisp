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

(defgeneric copy-genotypes (from to metadata &key key test)
  (:documentation "Copies genotype calls (AA, AB, BB) FROM sources
such as GTC files TO a destination such as a Plink BED file. The
METADATA argument is used to provide additional information required
during the operation.

Arguments:

- from (object): A producer of genotype data, such as a GTC object
  or a Lisp stream.
- to (object): A consumer of genotype data, such as a BED object or
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

(defmethod copy-genotypes ((gtc gtc) (bed bed) (snps vector) &key key test)
  (with-slots (stream)
      bed
    (let ((genotypes (data-field-of gtc :genotypes))
          (snps (if test
                    (remove-if test snps :key key)
                    snps)))
      (loop
         with selected = (make-array (length snps))
         for snp across snps
         for i = 0 then (1+ i)
         do (setf (aref selected i) (svref genotypes (1- (snp-index snp))))
         finally (write-bed-genotypes selected stream))))
  bed)

(defmethod copy-genotypes ((gtc gtc) (bed bed) (manifest bpm) &key key test)
  (copy-intensities gtc bed (snps-of manifest :key key :test test)))

(defgeneric gtc-to-bed (bed-filespec manifest gtc-filespecs &key test key)
  (:documentation "Creates a new Plink BED file containing the
aggregated genotype data from a list of GTC files. Also creates the
corresponding FAM and BIM annotation files.")
  (:method (bed-filespec (manifest bpm) gtc-filespecs &key test key)
    (with-bed (bed bed-filespec :direction :output :if-exists :supersede
                   :if-does-not-exist :create
                   :orientation :individual-major)
      (with-open-file (fam (plink-pathname bed-filespec "fam")
                           :direction :output :if-exists :supersede
                           :if-does-not-exist :create)
        (let ((snps (snps-of manifest :test test :key key))
              (prev-manifest-name))
          (with-open-file (bim (plink-pathname bed-filespec "bim")
                               :direction :output :if-exists :supersede
                               :if-does-not-exist :create)
            (loop
               for snp across snps
               do (write-bim-snp snp bim)))
          (dolist (gtc-filespec gtc-filespecs bed)
            (with-gtc (gtc gtc-filespec)
              (let ((manifest-name (data-field-of gtc :snp-manifest)))
                (if (null prev-manifest-name)
                    (setf prev-manifest-name manifest-name)
                    (check-arguments (string= prev-manifest-name manifest-name)
                                     (gtc-filespecs)
                                     "manifest ~s does not match previous ~
                                     manifest ~s" manifest-name
                                     prev-manifest-name)))
              (write-fam-individual (data-field-of gtc :sample-name) fam)
              (copy-genotypes gtc bed snps))))))))

(defun plink-pathname (bed-filespec type)
  "Returns the default Plink name for file of TYPE given a pathname
designator BED-FILESPEC."
  (check-arguments (string-equal "bed" (pathname-type bed-filespec))
                   (bed-filespec)
                   "not a BED filespec")
  (make-pathname :type type :defaults bed-filespec))

