;;;
;;; Copyright (c) 2012 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
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

(in-package :cl-user)

;; Assumes a recent ASDF with configuration API
(require 'asdf)

;; Set the output translation to put fasl files in the build directory
(asdf:initialize-output-translations
 (list :output-translations (list t (merge-pathnames "build/**/*.*"))
       :ignore-inherited-configuration))

(asdf:load-system :deoxybyte-systems)
(asdf:load-system :genotype-call)
(asdf:clear-configuration)

(sb-ext:save-lisp-and-die "genotype-call"
                          :executable t
                          :save-runtime-options t
                          :toplevel #'genotype-call:cli)
