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

(defparameter *commands* (make-hash-table :test #'equal)
  "A hash-table mapping a command name to a list containing a CLI
object and a function implementing the desired command.")

(defun get-cli (cmd)
  "Returns the CLI object for command string CMD."
  (first (gethash cmd *commands*)))

(defun get-fn (cmd)
  "Returns the function for command string CMD."
  (second (gethash cmd *commands*)))

(defun register-command (cmd cli fn)
  "Registers string CMD as calling FN via CLI. CLI must be a symbol
designating a CLI class."
  (check-arguments (and (stringp cmd) (symbolp cli) (functionp fn))
                   (cmd cli fn)
                   "expected a command string, a cli symbol and a function")
  (setf (gethash cmd *commands*) (list (make-instance cli) fn)))

(defun print-avail-commands (&optional (stream *error-output*))
  "Prints help for all registered commands to STREAM."
  (let ((commands (sort (loop
                           for cmd being the hash-keys of *commands*
                           collect cmd) #'string<)))
    (write-line "Available commands:" stream)
    (terpri stream)
    (mapc (lambda (cmd)
            (let* ((cli (get-cli cmd))
                   (doc (documentation-of cli)))
              (if doc
                  (help-message cli doc stream)
                  (warn "No help was found for ~a" cmd)))) commands)))

(define-cli input-mixin ()
  ((input "input" :required-option t :value-type 'string
          :documentation "The input file or stream.")))

(define-cli output-mixin ()
  ((output "output" :required-option t :value-type 'string
           :documentation "The output file or stream.")))

(define-cli manifest-mixin ()
  ((manifest "manifest" :required-option t :value-type 'string
             :documentation "The BeadPool Manifest file or stream.")))

(define-cli chromosome-mixin ()
  ((chromosome "chromosome" :required-option nil :value-type 'string
               :documentation "Limit processing to the named chromosome.")))

(define-cli gtc-to-sim-cli (cli output-mixin manifest-mixin chromosome-mixin)
  ()
  (:documentation "gtc-to-sim --output <filename> --manifest <filename>
[--chromosome <name>] <GTC filenames ...>"))

(define-cli sim-to-illuminus-cli (cli input-mixin output-mixin manifest-mixin
                                      chromosome-mixin)
  ((start "start" :required-option nil :value-type 'integer
          :documentation "The start of the range of SNPs to process.")
   (end "end" :required-option nil :value-type 'integer
        :documentation "The end of the range of SNPs to process."))
  (:documentation "sim-to-illuminus --input <filename> --output <filename>
--manifest <filename> [--chromosome <name>]"))

(define-cli mock-study-cli (cli manifest-mixin)
  ((study-name "study-name" :required-option t :value-type 'string
               :documentation "The name of the study, used in file naming.")
   (num-samples "num-samples" :required-option t :value-type 'integer
                :documentation "The number of samples and hence GTC files.")
   (num-snps "num-snps" :required-option t :value-type 'integer
             :documentation "The number of SNPs in the mock manifest."))
  (:documentation "mock-study --study-name <name> --num-samples <n>
--num-snps <m> --manifest <filename>"))

(defun cli ()
  "Applies the appropriate command line interface."
  (flet ((errmsg (c)
           (write-line (string-capitalize (format nil "~a" c) :end 1)
                       *error-output*)
           (terpri *error-output*)))
    (with-argv (argv)
      (let* ((cmd (first argv))
             (args (rest argv))
             (cli (get-cli cmd)))
        (handler-case
            (if (or (null cmd) (null cli))
                (error 'unknown-command :cli cli :command cmd)
                (multiple-value-bind (parsed unmatched)
                    (parse-command-line cli args)
                  (funcall (get-fn cmd) parsed unmatched)))
          (unknown-command (condition)
            (errmsg condition)
            (print-avail-commands)
            (quit-lisp :status 2))
          (cli-error (condition)
            (errmsg condition)
            (princ "Usage: " *error-output*)
            (let ((msg (documentation-of cli)))
              (if msg
                  (help-message cli msg *error-output*)
                  (warn "No help was found for ~a~%" cmd)))
            (quit-lisp :status 3))
          (file-error (condition)
            (errmsg condition))
          (error (condition)
            (errmsg condition)
            (write-line "Backtrace follows:" *error-output*)
            (error condition)))))))

(register-command "gtc-to-sim" 'gtc-to-sim-cli
                  (lambda (parsed-args &optional other)
                    (let ((output (option-value 'output parsed-args))
                          (manifest (load-bpm
                                     (option-value 'manifest parsed-args)))
                          (chromosome (option-value 'chromosome parsed-args)))
                      (if chromosome
                          (gtc-to-sim output manifest other
                                      :test (make-chromosome-p
                                             manifest chromosome #'string=)
                                      :key #'snp-chromosome)
                          (gtc-to-sim output manifest other)))))

(register-command "sim-to-illuminus" 'sim-to-illuminus-cli
                  (lambda (parsed-args &optional other)
                    (declare (ignorable other))
                    (let ((input (maybe-standard-stream
                                  (option-value 'input parsed-args)))
                          (output (maybe-standard-stream
                                   (option-value 'output parsed-args)))
                          (manifest (load-bpm
                                     (option-value 'manifest parsed-args)))
                          (chromosome (option-value 'chromosome parsed-args))
                          (start (or (option-value 'start parsed-args) 0))
                          (end (option-value 'end parsed-args)))
                      (if chromosome
                          (multiple-value-bind (cstart cend)
                              (chromosome-boundaries manifest chromosome)
                            (sim-to-illuminus output manifest input
                                              :start (+ cstart start)
                                              :end (min cend (+ cstart end))))
                          (sim-to-illuminus output manifest input
                                            :start start :end end)))))

(register-command "mock-study" 'mock-study-cli
                  (lambda (parsed-args &optional other)
                    (declare (ignorable other))
                    (let ((manifest (option-value 'manifest parsed-args)))
                      (generate-manifest
                       manifest (option-value 'num-snps parsed-args))
                      (generate-study
                       (option-value 'study-name parsed-args)
                       (option-value 'num-samples parsed-args)
                       manifest))))
