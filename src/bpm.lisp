;;;
;;; Copyright (c) 2011 Genome Research Ltd. All rights reserved.
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

(defvar *bpm-header*
  "Index,Name,Chromosome,Position,GenTrain Score,SNP,ILMN Strand,Customer Strand,NormID"
  "The header expected in a Beapool Manifest CSV file.")

(defparameter *default-bpm-size* 100000
  "The default size (number of SNPs) per manifest assumed prior to
  parsing. Vectors created while parsing start at this size.")

(defclass bpm ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The manifest input stream.")
   (name-index :initform (make-hash-table :test #'equal) :initarg :name-index
               :documentation "A mapping of SNP name to SNP index.")
   (chr-index :initform (make-hash-table :test #'equal) :initarg :chr-index
              :documentation "A mapping of chromosome name to a list
              of SNP indices for that chromosome.")
   (names :initform (make-array 0) :initarg :names :reader names-of
          :documentation "A vector of SNP names, in SNP order.")
   (positions :initform (make-array 0) :initarg :positions :reader positions-of
              :documentation "A vector of 1-based SNP chromosome
              positions, in SNP order.")
   (alleles :initform (make-array 0) :initarg :alleles :reader alleles-of)
   (ilmn-strand :initform (make-array 0) :initarg :ilmn-strand
                :reader ilmn-strand-of)
   (cust-strand :initform (make-array 0) :initarg :cust-strand
                :reader cust-strand-of)
   (normalization :initform (make-array 0) :initarg :normalization
                  :reader normalization-of))
  (:documentation "An Illumina Beadpool Manifest. For convenience,
  SNPs are indexed here starting from 0, rather than from 1 as in the
  manifest file."))

(defmethod print-object ((manifest bpm) stream)
  (print-unreadable-object (manifest stream :type t :identity nil)
    (with-accessors ((s stream-of) (names names-of))
        manifest
      (format stream "~@[~a, ~]~d SNPs"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s)) (length names)))))

(defgeneric chromosomes-of (manifest)
  (:documentation "Returns a sorted list of the chromosomes described
  by MANIFEST.")
  (:method ((manifest bpm))
    (with-slots (chr-index)
        manifest
      (loop
         for chr being the hash-keys of chr-index
         collect chr into chrs
         finally (return (sort chrs #'string<))))))

(defgeneric num-snps-of (manifest &optional chromosome)
  (:documentation "Returns the number of SNPs described by MANIFEST,
  optionally restricting the count to only those SNPs on CHROMOSOME.")
  (:method ((manifest bpm) &optional chromosome)
    (with-slots (name-index chr-index)
        manifest
      (if chromosome
          (with-accessors ((chromosomes chromosomes-of))
              manifest
            (check-arguments (member chromosome chromosomes :test #'string=)
                             (chromosome)
                             "invalid chromosome, expected one of ~a"
                             chromosomes)
            (length (gethash chromosome chr-index)))
          (hash-table-count name-index)))))

(defgeneric find-snp (manifest identifier)
  (:documentation "Returns a SNP name given a SNP index as IDENTIFIER,
  or a SNP index given a SNP name as IDENTIFIER.")
  (:method ((manifest bpm) (name string))
    (gethash name (slot-value manifest 'name-index)))
  (:method ((manifest bpm) (index fixnum))
    (aref (slot-value manifest 'names) index)))

(defun read-bpm (stream)
  "Returns a new BPM object from read STREAM."
  (check-arguments (streamp stream) (stream) "expected a stream argument")
  (let ((header (read-line stream t :eof)))
    (unless (starts-with-string-p header *bpm-header*)
      (error 'malformed-file-error
             :format-control "invalid header: expected ~s, found ~s"
             :format-arguments (list *bpm-header* header)))
    (let ((name-index (make-hash-table :test #'equal :size *default-bpm-size*))
          (chr-index (make-hash-table :test #'equal))
          (names (make-array *default-bpm-size* :element-type 'string
                             :adjustable t :fill-pointer 0))
          (posn (make-array *default-bpm-size* :element-type 'fixnum
                            :adjustable t :fill-pointer 0))
          (norm (make-array *default-bpm-size* :element-type 'fixnum
                            :adjustable t :fill-pointer 0)))
      (flet ((index-snp (rec i)
               (setf (gethash (assocdr 'name rec) name-index) i))
             (index-chr (rec i)
               (let ((chr (assocdr 'chromosome rec)))
                 (push i (gethash chr chr-index (list)))))
             (store (rec key vector)
               (vector-push-extend (assocdr key rec) vector) 1000)
             (make-bpm (i)
               (make-instance
                'bpm :stream stream :name-index name-index :chr-index chr-index
                :names (make-array i :element-type 'string
                                   :initial-contents names)
                :positions (make-array i :element-type 'fixnum
                                       :initial-contents posn)
                :normalization (rank-norm-ids
                                (make-array i :element-type 'fixnum
                                            :initial-contents norm)))))
        (do ((line (read-line stream nil :eof) (read-line stream nil :eof))
             (i 0 (1+ i)))
            ((eql :eof line) (make-bpm i))
          (let ((record (parse-bpm-record line)))
            (unless (= (1+ i) (assocdr 'index record))
              (error 'malformed-file-error
                     "SNP records were not in ascending order"))
            (index-snp record i)
            (index-chr record i)
            (store record 'name names)
            (store record 'position posn)
            (store record 'norm-id norm)))))))

(define-line-parser parse-bpm-record #\,
  ((index :type :integer)
   (name :type :string)
   (chromosome :type :string)
   (position :type :integer)
   (gentrain-score :type :float)
   (snp :type :string)
   (ilmn-strand :type :string :parser #'parse-strand)
   (customer-strand :type :string :parser #'parse-strand)
   (norm-id :type :integer)))

(declaim (inline parse-strand))
(defun parse-strand (field-name str &key (start 0) end null-str)
  (declare (ignore field-name null-str))
  (declare (optimize (speed 3)))
  (declare (type simple-string str))
  (intern (subseq str start end) :keyword))

(defun rank-norm-ids (norm-ids)
  "Modifies vector NORM-IDS containing all NormIDs in the manifest,
by substituting each NormID with its rank on the sorted set of all
unique NormIDs. Ranks start from 0. The rank for each SNP is then an
index into the vector of XForms within any GTC file and indicates
which normalization parameters (XForm) is to be used for each SNP."
  (let* ((rank-order (sort (remove-duplicates norm-ids) #'<))
         (rank-map (make-hash-table :size (length rank-order))))
    (loop
       for rank from 0 below (length rank-order)
       do (setf (gethash (aref rank-order rank) rank-map) rank))
    (loop
       for i from 0 below (length norm-ids)
       do (setf (aref norm-ids i) (gethash (aref norm-ids i) rank-map))
       finally (return norm-ids))))
