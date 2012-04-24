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

(defvar *gtc-magic* (map-into (make-array 3 :element-type 'octet)
                              #'char-code "gtc")
  "Genotype Call (GTC) file magic string.")

(defvar *gtc-version* 3
  "Genotype Call (GTC) file version.")

(defclass gtc ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The input or output stream.")
   (version :initform nil :reader version-of
            :documentation "The GTC format version.")
   (toc :initform nil :reader toc-of
        :documentation "The table of contents. A vector of TOC-ENTRY
        objects.")
   (manifest-name :initform nil :reader manifest-name-of
                  :documentation "A cached value of the SNP manifest
                  name for internal use."))
  (:documentation "An Illumina Genotype Call binary file. These files
  are created by the Illumina genotyping process."))

(defclass toc-entry ()
  ((name :initform :unknown :initarg :name :reader name-of
         :documentation "The data name, a symbol.")
   (position :initarg :position :reader position-of
             :documentation "The position in the stream at which the
             data lie.")
   (parser :initform nil :initarg :parser :reader parser-of
           :documentation "The parser function for the data.")
   (immediate :initform nil :initarg :immediate :reader immediate-value-p))
  (:documentation "An Illumina Genotype Call file table of contents
  entry."))

(defstruct (xform (:conc-name xf-)
                  (:constructor make-xform (&optional
                                            version offset-x offset-y
                                            scale-x scale-y
                                            shear theta
                                            reserved0 reserved1
                                            reserved2 reserved3
                                            reserved4 reserved5)))
  (version 0 :type uint32 :read-only t)
  (offset-x 0.f0 :type single-float :read-only t)
  (offset-y 0.f0 :type single-float :read-only t)
  (scale-x 0.f0 :type single-float :read-only t)
  (scale-y 0.f0 :type single-float :read-only t)
  (shear 0.f0 :type single-float :read-only t)
  (theta 0.f0 :type single-float :read-only t)
  (reserved0 0.f0 :type single-float :read-only t)
  (reserved1 0.f0 :type single-float :read-only t)
  (reserved2 0.f0 :type single-float :read-only t)
  (reserved3 0.f0 :type single-float :read-only t)
  (reserved4 0.f0 :type single-float :read-only t)
  (reserved5 0.f0 :type single-float :read-only t))

(defmethod initialize-instance :after ((gtc gtc) &key)
  (with-slots (stream version toc manifest-name)
      gtc
    (when (input-stream-p stream)
      (let ((buffer (make-array 4 :element-type 'octet :initial-element 0)))
        (read-magic stream buffer *gtc-magic*)
        (setf version (read-version stream buffer)
              toc (read-gtc-toc stream buffer)
              manifest-name (data-field-of gtc :snp-manifest))))))

(defmethod print-object ((gtc gtc) stream)
  (print-unreadable-object (gtc stream :type t :identity nil)
    (with-accessors ((s stream-of) (version version-of) (toc toc-of))
        gtc
      (format stream "~@[~a ~][GTC version ~d] contents:~%  ~a~%"
              (when (subtypep (type-of s) 'file-stream)
                (file-namestring s)) version toc))))

(defmethod print-object ((entry toc-entry) stream)
  (print-unreadable-object (entry stream :type t :identity nil)
    (format stream "~a ~d" (name-of entry) (position-of entry))))

(defgeneric data-field-of (gtc name)
  (:documentation "Reads data from GTC and returns a parsed value. The NAME
argument is a keyword symbol indicating which data to return.

Arguments:

- gtc (object): A gtc instance.
- name (symbol): A symbol, one of :num-snps :sample-name :sample-plate
  :sample-well :cluster-file :snp-manifest:imaging-date :autocall-date
  :autocall-version :normalization-xforms :x-controls :y-controls
  :x-intensities :y-intensities :genotypes :basecalls :genotype-scores
  :scanner-data

Returns:

- Data value (type depends on data)")
  (:method ((gtc gtc) (name symbol))
    (let ((toc-entry (find name (toc-of gtc) :key #'name-of)))
      (check-arguments toc-entry (name) "data are not present in this GTC file")
      (read-gtc-data-field gtc toc-entry))))

(defmacro with-gtc ((var filespec &rest args) &body body)
  `(let ((,var (gtc-open ,filespec ,@args)))
     (unwind-protect
          (progn
            ,@body)
       (when ,var
         (gtc-close ,var)))))

(defun gtc-open (filespec &rest args)
  (let ((stream (apply #'open filespec :element-type 'octet args)))
    (if (open-stream-p stream)
        (make-instance 'gtc :stream stream)
        (error 'invalid-operation-error
               :format-control "GTC stream ~a closed unexpectedly"
               :format-arguments (list stream)))))

(defun gtc-close (gtc)
  (close (stream-of gtc)))

(defun normalize (x y xform)
  "Returns X and Y intensities as two values after normalization with
alist XFORM."
  (declare (optimize (speed 3)))
  (let ((offset-x (xf-offset-x xform))
        (offset-y (xf-offset-y xform))
        (scale-x (xf-scale-x xform))
        (scale-y (xf-scale-y xform))
        (shear (xf-shear xform))
        (theta (xf-theta xform)))
    (declare (type uint16 x y)
             (type single-float offset-x offset-y scale-x scale-y shear theta))
    (let* ((x1 (- x offset-x))
           (y1 (- y offset-y))
           (x2 (+ (* (cos theta) x1) (* (sin theta) y1)))
           (y2 (+ (* (- (sin theta)) x1) (* (cos theta) y1)))
           (x3 (- x2 (* shear y2)))
           (normx (/ x3 scale-x))
           (normy (/ y2 scale-y)))
      ;; Eliminate negative values
      (values (if (minusp normx) 0.f0 normx)
              (if (minusp normy) 0.f0 normy)))))

(defun read-gtc-toc (stream buffer)
  "Reads the GTC table of contents fom STREAM as a simple-array."
  (let ((toc (make-array (read-uint32 stream buffer) :initial-element nil)))
    (loop
       for i from 0 below (length toc)
       do (setf (aref toc i) (read-gtc-toc-entry stream buffer))
       finally (return toc))))

(defun read-gtc-toc-entry (stream buffer)
  "Reads a single toc-entry from STREAM."
  (let* ((id (read-uint16 stream buffer))
         (position (read-uint32 stream buffer)))
    (destructuring-bind (name fname &optional immediate)
        (case id
          (1 '(:num-snps read-uint32 t))
          (10 '(:sample-name read-gtc-string))
          (11 '(:sample-plate read-gtc-string))
          (12 '(:sample-well read-gtc-string))
          (100 '(:cluster-file read-gtc-string))
          (101 '(:snp-manifest read-gtc-string))
          (200 '(:imaging-date read-gtc-string))
          (201 '(:autocall-date read-gtc-string))
          (300 '(:autocall-version read-gtc-string))
          (400 '(:normalization-xforms read-gtc-xforms))
          (500 '(:x-controls read-gtc-intensities))
          (501 '(:y-controls read-gtc-intensities))
          (1000 '(:x-intensities read-gtc-intensities))
          (1001 '(:y-intensities read-gtc-intensities))
          (1002 '(:genotypes read-gtc-genotypes))
          (1003 '(:basecalls read-gtc-basecalls))
          (1004 '(:genotype-scores read-gtc-genotype-scores))
          (1005 '(:scanner-data read-gtc-scanner-data))
          (otherwise '(:unknown nil)))
      (make-instance 'toc-entry :name name
                     :parser (and fname (symbol-function fname))
                     :position position :immediate immediate))))

(defun write-gtc-toc-entry (name value stream buffer)
  "Writes a single toc entry identified by symbol NAME and having
uint32 VALUE to STREAM (see GTC specification for valid TOC names and
values)."
  (let ((id (ecase name
              (:num-snps 1)
              (:sample-name 10)
              (:snp-manifest 101)
              (:normalization-xforms 400)
              (:x-intensities 1000)
              (:y-intensities 1001)
              (:genotypes 1002)
              (:basecalls 1003)
              (:genotype-scores 1004))))
    (write-uint16 id stream buffer)
    (write-uint32 value stream buffer)
    name))

(defun read-gtc-data-field (gtc toc-entry)
  "Returns a parsed data element denoted by TOC-ENTRY in the GTC toc."
  (with-slots (stream)
      gtc
    (if (immediate-value-p toc-entry)
        (position-of toc-entry)
        (let ((buffer (make-array 4 :element-type 'octet :initial-element 0))
              (fn (parser-of toc-entry)))
          (with-restored-position stream
            (file-position stream (position-of toc-entry))
            (funcall fn stream buffer))))))

(defun read-gtc-intensities (stream buffer)
  "Reads intensities from STREAM as a vector of uint16."
  (declare (optimize (speed 3)))
  (let* ((n (read-uint32 stream buffer))
         (num-bytes (* 2 (the uint32 n)))
         (intensities (make-array n :element-type 'uint16 :initial-element 0))
         (read-buffer (read-record
                       stream (make-array num-bytes :element-type 'octet
                                          :initial-element 0) num-bytes)))
    (loop
       for i from 0 below n
       for j from 0 below num-bytes by 2
       do  (setf (aref intensities i) (decode-uint16le read-buffer j))
       finally (return intensities))))

(defun write-gtc-intensities (intensities stream buffer)
  "Writes vector of uint16 INTENSITIES to STREAM."
  (write-uint32 (length intensities) stream buffer)
  (loop
     for x across intensities
     do (write-uint16 x stream buffer))
  intensities)

(defun read-gtc-xforms (stream buffer)
  "Reads intensity normalization transforms from STREAM as a
simple-array of alists. Not all of the transforms are applied to every
SNP intensity; the information on which to apply is contained in the
beadpool manifest."
  (let* ((n (read-uint32 stream buffer))
         (xforms (make-array n :initial-element nil)))
    (loop
       for i from 0 below n
       do (let ((version (read-uint32 stream buffer))
                (fields (loop
                           repeat 12
                           collect (read-float stream buffer))))
            (setf (aref xforms i) (apply #'make-xform version fields)))
       finally (return xforms))))

(defun write-gtc-xforms (xforms stream buffer)
  (write-uint32 (length xforms) stream buffer)
  (loop
     for xform across xforms
     do (progn
          (write-uint32 (xf-version xform) stream buffer)
          (dolist (reader '(xf-offset-x xf-offset-y
                            xf-scale-x xf-scale-y
                            xf-shear xf-theta
                            xf-reserved0 xf-reserved1
                            xf-reserved2 xf-reserved3
                            xf-reserved4 xf-reserved5))
            (write-float (funcall reader xform) stream buffer))))
  xforms)

(defun read-gtc-genotypes (stream buffer)
  "Reads genotype calls from STREAM as a simple-vector of strings."
  (let* ((len (read-uint32 stream buffer))
         (genotypes (make-array len :element-type t :initial-element nil))
         (bytes (read-record stream (make-array len :element-type 'octet
                                                :initial-element 0) len)))
    (loop
       for i from 0 below len
       do (setf (svref genotypes i) (ecase (aref bytes i)
                                      (0 nil)
                                      (1 "AA")
                                      (2 "AB")
                                      (3 "BB")))
       finally (return genotypes))))

(defun write-gtc-genotypes (genotypes stream buffer)
  "Writes simple-vector of strings GENOTYPES to STREAM."
  (write-uint32 (length genotypes) stream buffer)
  (let ((bytes (make-array (length genotypes) :element-type 'octet
                           :initial-element 0)))
    (loop
       for i from 0 below (length genotypes)
       do (setf (aref bytes i) (let ((g (svref genotypes i)))
                                 (cond ((string= "AA" g)
                                        1)
                                       ((string= "AB" g)
                                        2)
                                       ((string= "BB" g)
                                        3)
                                       (t
                                        0)))))
    (write-sequence bytes stream)))

(defun read-gtc-genotype-scores (stream buffer)
  "Reads genotype scores from STREAM as a vector of single-floats."
  (let* ((n (read-uint32 stream buffer))
         (scores (make-array n :element-type 'float :initial-element 0.f0)))
    (loop
         for i from 0 below n
         do (setf (aref scores i) (read-float stream buffer))
         finally (return scores))))

(defun write-gtc-genotype-scores (scores stream buffer)
  "Writes vector of single-floats SCORES to STREAM."
  (write-uint32 (length scores) stream buffer)
  (loop
     for score across scores
     do (write-float score stream buffer))
  scores)

(defun read-gtc-basecalls (stream buffer)
  "Reads basecalls from STREAM as a string."
  (let* ((n (read-uint32 stream buffer))
         (basecalls (make-array n :initial-element nil)))
    (loop
       for i below n
       do (let ((rec (read-record stream buffer 2)))
            (setf (svref basecalls i) (octets-to-string rec 0 2)))
       finally (return basecalls))))

(defun write-gtc-basecalls (basecalls stream buffer)
  "Writes simple-vector of strings BASECALLS to STREAM."
  (write-uint32 (length basecalls) stream buffer)
  (let ((bytes (make-array (* 2 (length basecalls)) :element-type 'octet)))
    (loop
       for i from 0 below (length basecalls)
       for j from 0 below (length bytes) by 2
       do (let ((calls (svref basecalls i)))
            (setf (aref bytes j) (char-code (char calls 0))
                  (aref bytes (1+ j)) (char-code (char calls 1)))))
    (write-sequence bytes stream)))

(defun read-gtc-scanner-data (stream buffer)
  "Reads the scanner data from STREAM as an alist."
  (let* ((scanner-name (read-gtc-string stream buffer))
         (pmt-green (read-uint32 stream buffer))
         (pmt-red (read-uint32 stream buffer))
         (scanner-version (read-gtc-string stream buffer))
         (imaging-user (read-gtc-string stream buffer)))
    (pairlis '(:scanner-name :pmt-green :pmt-red
               :scanner-version :imaging-user)
             (list scanner-name pmt-green pmt-red
                   scanner-version imaging-user))))
