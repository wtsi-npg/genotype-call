;;;
;;; Copyright (c) 2011-2012 Genome Research Ltd. All rights reserved.
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
  "The header expected in a Beadpool Manifest CSV file.")

(defparameter *default-bpm-size* 100000
  "The default size (number of SNPs) per manifest assumed prior to
  parsing. Vectors created while parsing start at this size.")

(deftype snp-count ()
  "SNP index type denoting the possible range of number of SNPs per
  microarray. This value is set at an order of magnitiude larger than
  the largest chip currently available."
  '(and fixnum (integer 0 100000000)))

(defstruct (snp (:constructor make-snp (index name chromosome position
                                        alleles ilmn-strand cust-strand
                                        norm-id)))
  "A representation of the data for a single SNP as stored in a
Beadpool Manifest."
  (index 0 :type fixnum)
  (name "" :type simple-string)
  (chromosome "" :type simple-string)
  (position 0 :type fixnum)
  (alleles "" :type simple-string)
  (ilmn-strand #\? :type character)
  (cust-strand #\? :type character)
  (norm-id 0 :type fixnum)
  (norm-rank 0 :type fixnum))

(defclass bpm ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The manifest input stream.")
   (snps :initform (make-array 0) :initarg :snps)
   (chromosomes :initform nil :initarg :chromosomes :reader chromosomes-of))
  (:documentation "An Illumina Beadpool Manifest."))

(defmethod initialize-instance :after ((manifest bpm) &key)
  (with-slots (snps)
      manifest
    (setf snps (stable-sort
                (normalize-snp-alleles (rank-norm-ids snps)) #'location<))))

(defmethod print-object ((manifest bpm) stream)
  (print-unreadable-object (manifest stream :type t :identity nil)
    (with-slots ((s stream) snps)
        manifest
      (format stream "~@[~a, ~]~d SNPs"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s)) (length snps)))))

(defgeneric snps-of (manifest &key key test)
  (:documentation "Returns the SNPs described by MANIFEST, optionally
  restricting the count to only those SNPs for which predicate TEST
  returns T.")
  (:method ((manifest bpm) &key key test)
    (with-slots (snps)
        manifest
      (if (null test)
          snps
          (let ((key (cond ((null key)
                            #'identity)
                           ((functionp key)
                            key)
                           (t
                            (fdefinition key))))
                (test (cond ((functionp test)
                             test)
                            (t
                             (fdefinition test)))))
            (loop
               for snp across snps
               when (funcall test (funcall key snp))
               count snp into n
               and collect snp into subset
               finally (return (make-array n :initial-contents subset))))))))

(defgeneric num-snps-of (manifest &key key test)
  (:documentation "Returns the number of SNPs described by MANIFEST,
  optionally restricting the count to only those SNPs for which
  predicate TEST returns T.")
  (:method ((manifest bpm) &key key test)
    (length (snps-of manifest :key key :test test))))

(defun cnv-probe-p (probe-name)
  "Returns T if the PROBE-NAME indicates a copy-number variation (CNV)
probe."
  (starts-with-string-p probe-name "cnv"))

(defun make-chromosome-p (manifest chromosome predicate)
  "Makes a new predicate that accepts a CHROMOSOME as an argument and
returns the result of testing CHROMOSOME with PREDICATE. CHROMOSOME
must be represneted in MANIFEST."
  (with-slots (chromosomes)
      manifest
    (check-arguments (member chromosome chromosomes :test predicate)
                     (chromosome)
                     "invalid chromosome, expected one of ~a"
                     chromosomes))
  (let ((predicate (if (functionp predicate)
                       predicate
                       (fdefinition predicate))))
    (lambda (arg)
      (funcall predicate chromosome arg))))

(defun load-bpm (filespec &key strict-ordering)
  (with-open-file (stream filespec :external-format :ascii)
    (read-bpm stream :strict-ordering strict-ordering)))

(defun read-bpm (stream &key (strict-ordering t))
  "Returns a new BPM object from read STREAM."
  (check-arguments (streamp stream) (stream) "expected a stream argument")
  (let ((header (read-line stream t :eof)))
    (unless (starts-with-string-p header *bpm-header*)
      (error 'malformed-file-error
             :format-control "invalid header: expected ~s, found ~s"
             :format-arguments (list *bpm-header* header)))
    (let ((snps (make-array *default-bpm-size* :adjustable t :fill-pointer 0))
          (chromosomes (make-hash-table :test #'equal)))
      (flet ((make-bpm ()
               (make-instance
                'bpm :stream stream
                :snps (make-array (length snps) :initial-contents snps)
                :chromosomes (loop
                                for chr being the hash-keys of chromosomes
                                collect chr into chrs
                                finally (return (sort chrs #'string<))))))
        (do ((line (read-line stream nil :eof) (read-line stream nil :eof))
             (i 0 (1+ i)))
            ((eql :eof line) (make-bpm))
          (let* ((record (parse-bpm-record line))
                 (index (assocdr 'index record))
                 (chromosome (assocdr 'chromosome record)))
            (when (and strict-ordering (/= (1+ i) index))
              (error 'malformed-file-error :file stream
                     :format-control "SNP records were not in ascending order ~
                                      with strict ordering, expected record ~
                                      ~d, found ~d"
                     :format-arguments (list (1+ i) index)))
            (setf (gethash chromosome chromosomes) t)
            (vector-push-extend
             (make-snp index
                       (assocdr 'name record)
                       chromosome
                       (assocdr 'position record)
                       (assocdr 'alleles record)
                       (assocdr 'ilmn-strand record)
                       (assocdr 'cust-strand record)
                       (assocdr 'norm-id record)) snps 50000)))))))

(define-line-parser parse-bpm-record #\,
  ((index :type :integer)
   (name :type :string)
   (chromosome :type :string)
   (position :type :integer)
   (gentrain-score :type :float :ignore t)
   (alleles :type :string :parser #'parse-alleles :validator #'valid-alleles-p)
   (ilmn-strand :type :string :parser #'parse-strand)
   (cust-strand :type :string :parser #'parse-strand)
   (norm-id :type :integer)))

(declaim (inline parse-strand))
(defun parse-strand (field-name str &key (start 0) end null-str)
  "Returns a character representing the SNP strand, one of T (TOP),
B (BOT), - (MINUS) or + (PLUS)."
  (declare (ignore field-name null-str))
  (declare (optimize (speed 3)))
  (declare (type simple-string str))
  (cond ((char= #\T (char str start))       ; TOP or Top
         #\T)
        ((char= #\B (char str start))       ; BOT
         #\B)
        ((char= #\P (char str start))       ; PLUS or P
         #\+)
        ((char= #\M (char str start))       ; MINUS or M
         #\-)
        (t
         (error 'malformed-field-error :field (subseq str start end)
                :record str))))

(defun parse-alleles (field-name str &key (start 0) end null-str)
  (declare (ignore field-name null-str))
  (cond ((string= "[N/A]" str :start2 start :end2 end)
         "NA")
        ((>= (length str) (+ start 3))
         (coerce (list (char str (1+ start)) (char str (+ start 3)))
                 'simple-string))
        (t
         (error 'malformed-field-error :field (subseq str start end)
                :record str))))

(defun nucleotidep (c)
  "Returns T if character C represents a nucleotide base."
  (or (char= #\A c) (char= #\C c) (char= #\G c) (char= #\T c)))

(defun complement-allele (c)
  "Returns a character representing the allele on the complementary
strand to allele character C."
  (ecase c
    (#\A #\T)
    (#\C #\G)
    (#\G #\C)
    (#\T #\A)
    (#\D #\I)
    (#\I #\D)))

(defun valid-alleles-p (str)
  "Returns T if STR is a valid manifest allele string."
  (or (string= "NA" str) (string= "DI" str) (string= "ID" str)
      (every #'nucleotidep str)))

;; The complete set of Ilmn_strand values across all the manifests that I
;; have available is somewhat variable, being:
;; { Bot BOT M MINUS P PLUS Top TOP }
;;
;; The complete set of allele values across all the manifests that I
;; have available is:
;; { [A/A] [A/C] [A/G] [A/T] [C/C] [C/G] [G/C] [G/G]
;;   [T/A] [T/C] [T/G] [T/T]
;;   [D/I] [I/D] [N/A] }

(defun normalize-snp-alleles (snps)
  "Modifies SNPs by applying NORMALIZE-ALLELES to each one."
  (loop
     for snp across snps
     do (let ((norm (normalize-alleles (snp-alleles snp)
                                       (snp-ilmn-strand snp))))
          (if (string= (snp-alleles snp) norm)
              snp
              (setf (snp-alleles snp) norm
                    (snp-ilmn-strand snp) #\T)))
     finally (return snps)))

(defun normalize-alleles (alleles strand)
  "Returns a string representing ALLELES, normalizing to the TOP
strand."
  (if (string= "NA" alleles)
      alleles
      (ecase strand
        ((#\T #\+ #\-) alleles)
        (#\B (map 'simple-string #'complement-allele alleles)))))

(defun location< (snp1 snp2)
  (let ((chr1 (snp-chromosome snp1))
        (chr2 (snp-chromosome snp2)))
    (if (string= chr1 chr2)
        (< (snp-position snp1) (snp-position snp2))
        (string< chr1 chr2))))

(defun rank-norm-ids (snps)
  "Modifies vector SNPS containing all SNPs in the manifest,
by setting the NORM-RANK field to a value equal to the rank of the
SNP's NORM-ID on the sorted set of all unique NORM-IDs. Ranks start
from 0. The rank for each SNP is then an index into the vector of
XForms within any GTC file and indicates which normalization
parameters (XForm) is to be used for each SNP."
  (let* ((norm-ids (loop
                      with ids = (make-hash-table)
                      for snp across snps
                      do (setf (gethash (snp-norm-id snp) ids) t)
                      finally (return (loop
                                         for id being the hash-keys of ids
                                         collect id))))
         (rank-order (sort (make-array (length norm-ids) :element-type 'fixnum
                                       :initial-contents norm-ids) #'<))
         (rank-map (loop
                      with ranks = (make-hash-table)
                      for rank from 0 below (length rank-order)
                      do (setf (gethash (aref rank-order rank) ranks) rank)
                      finally (return ranks))))
    (loop
       for rank from 0 below (length rank-order)
       do (setf (gethash (aref rank-order rank) rank-map) rank))
    (loop
       for snp across snps
       do (setf (snp-norm-rank snp) (gethash (snp-norm-id snp) rank-map))
       finally (return snps))))
