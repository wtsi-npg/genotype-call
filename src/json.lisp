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

(defmacro with-underscore-translation (&body body)
  "Sets up identifier transltion between Lisp and JSON."
  `(let ((json:*lisp-identifier-name-to-json* 'lisp-to-underscore)
         (json:*json-identifier-name-to-lisp* 'underscore-to-lisp))
     ,@body))

(defun lisp-to-underscore (str)
  (string-downcase
   (if (> (length str) 2)
       (substitute-if #\_ (lambda (c)
                            (char= #\- c)) str :start 1 :end (1- (length str)))
       str)))

(defun underscore-to-lisp (str)
  (string-upcase
   (if (> (length str) 2)
       (substitute-if #\- (lambda (c)
                            (char= #\_ c)) str :start 1 :end (1- (length str)))
       str)))

(defun read-json-sample-specs (stream)
  (with-underscore-translation
    (mapcar #'ensure-sample-urn (json:decode-json stream))))

(defun ensure-sample-urn (spec)
  (let ((uri (assocdr :uri spec))
        (spec-str (json:encode-json-to-string spec)))
    (check-field uri spec-str uri
                 "the sample URI is missing from this record")
    (let ((urn (puri:uri uri)))
      (check-field (urn-p urn) spec-str uri
                   "the URN in this record is malformed")
      (subst urn uri spec :test #'equal))))
