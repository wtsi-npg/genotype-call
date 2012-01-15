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

(defvar *gtc-magic* (map-into (make-array 3 :element-type 'octet)
                              #'char-code "gtc")
  "Genotype Call (GTC) file magic string.")

(defvar *xform-fields* '(:offset-x :offset-y
                         :scale-x :scale-y :shear :theta
                         :reserved0 :reserved1 :reserved2
                         :reserved3 :reserved4 :reserved5)
  "The fields of a normalization XForm.")

(defclass gtc ()
  ((stream :initform nil :initarg :stream :reader stream-of
           :documentation "The input stream.")
   (version :initform nil :reader version-of
            :documentation "The GTC format version.")
   (toc :initform nil :reader toc-of
        :documentation "The table of contents. A vector of TOC-ENTRY
        objects."))
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

(defmethod initialize-instance :after ((gtc gtc) &key)
  (with-slots (stream version toc)
      gtc
    (let ((buffer (make-array 4 :element-type 'octet :initial-element 0)))
      (read-magic stream buffer *gtc-magic*)
      (setf version (read-version stream buffer)
            toc (read-gtc-toc stream buffer)))))

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
      (check-arguments toc-entry (name) "is not a GTC data field")
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
               :format-control "GTC stream ~a  closed unexpectedly"
               :format-arguments (list stream)))))

(defun gtc-close (gtc)
  (close (stream-of gtc)))

(defun normalize (x y xform)
  "Returns X and Y intensities as two values after normalization with
alist XFORM."
  (declare (optimize (speed 3)))
  (let ((offset-x (assocdr :offset-x xform))
        (offset-y (assocdr :offset-y xform))
        (scale-x (assocdr :scale-x xform))
        (scale-y (assocdr :scale-y xform))
        (shear (assocdr :shear xform))
        (theta (assocdr :theta xform)))
    (declare (type uint16 x y)
             (type single-float offset-x offset-y scale-x scale-y shear theta))
    (let* ((x1 (- x offset-x))
           (y1 (- y offset-y))
           (x2 (+ (* (cos theta) x1) (* (sin theta) y1)))
           (y2 (+ (* (- (sin theta)) x1) (* (cos theta) y1)))
           (x3 (- x2 (* shear y2))))
      (values (/ x3 scale-x) (/ y2 scale-y)))))

(defun read-gtc-toc (stream buffer)
  "Reads the GTC table of contents fom STREAM as a simple-array."
  (let ((toc (make-array (decode-uint32le (read-record stream buffer 4))
                         :initial-element nil)))
    (loop
       for i from 0 below (length toc)
       do (setf (aref toc i) (read-gtc-toc-entry stream buffer))
       finally (return toc))))

(defun read-gtc-toc-entry (stream buffer)
  "Reads a single toc-entry from STREAM."
  (let* ((id (decode-uint16le (read-record stream buffer 2)))
         (position (decode-uint32le (read-record stream buffer 4))))
    (destructuring-bind (name fname &optional immediate)
        (case id
          (1 '(:num-snps read-uint32 t))
          (10 '(:sample-name read-string))
          (11 '(:sample-plate read-string))
          (12 '(:sample-well read-string))
          (100 '(:cluster-file read-string))
          (101 '(:snp-manifest read-string))
          (200 '(:imaging-date read-string))
          (201 '(:autocall-date read-string))
          (300 '(:autocall-version read-string))
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
         (intensities (make-array n :element-type 'uint16
                                  :initial-element 0))
         (read-buffer (read-record
                       stream (make-array num-bytes :element-type 'octet
                                          :initial-element 0) num-bytes)))
    (loop
       for i from 0 below n
       for j from 0 below num-bytes by 2
       do  (setf (aref intensities i) (decode-uint16le read-buffer j))
       finally (return intensities))))

(defun read-gtc-xforms (stream buffer)
  "Reads intensity normalization transforms from STREAM as a
simple-array of alists. Not all of the transforms are applied to every
SNP intensity; the information on which to apply is contained in the
beadpool manifest."
  (let* ((n (read-uint32 stream buffer))
         (xforms (make-array n :initial-element nil)))
    (loop
       for i from 0 below n
       do (let ((version (read-uint32 stream buffer)))
            (declare (ignore version))
            (setf (aref xforms i)
                  (mapcar (lambda (key)
                            (cons key (read-float stream buffer)))
                          *xform-fields*)))
         finally (return xforms))))

(defun read-gtc-genotypes (stream buffer)
  "Reads genotype calls from STREAM as a simple-array of strings."
  (let* ((len (read-uint32 stream buffer))
         (genotypes (make-array len :element-type t :initial-element nil))
         (bytes (read-record stream  (make-array len :element-type 'octet
                                                 :initial-element 0) len)))
    (loop
       for i from 0 below len
       do (setf (aref genotypes i) (ecase (aref bytes i)
                                     (0 nil)
                                     (1 "AA")
                                     (2 "AB")
                                     (3 "BB")))
       finally (return genotypes))))

(defun read-gtc-genotype-scores (stream buffer)
  "Reads genotype scores from STREAM as a vector of single-floats."
  (let* ((n (read-uint32 stream buffer))
         (scores (make-array n :element-type 'float :initial-element 0.f0)))
    (loop
         for i from 0 below n
         do (setf (aref scores i) (read-float stream buffer))
         finally (return scores))))

(defun read-gtc-basecalls (stream buffer)
  "Reads basecalls from STREAM as a string."
  (let* ((n (read-uint32 stream buffer))
         (basecalls (make-array n :initial-element nil)))
    (loop
       for i below n
       do (let ((rec (read-record stream buffer 2)))
            (setf (aref basecalls i) (octets-to-string rec 0 2)))
       finally (return basecalls))))

(defun read-gtc-scanner-data (stream buffer)
  "Reads the scanner data from STREAM as an alist."
  (let* ((scanner-name (read-string stream buffer))
         (pmt-green (read-uint32 stream buffer))
         (pmt-red (read-uint32 stream buffer))
         (scanner-version (read-string stream buffer))
         (imaging-user (read-string stream buffer)))
    (pairlis '(:scanner-name :pmt-green :pmt-red
               :scanner-version :imaging-user)
             (list scanner-name pmt-green pmt-red
                   scanner-version imaging-user))))
