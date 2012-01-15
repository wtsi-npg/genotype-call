;;;
;;; Copyright (c) 2011 Genome Research Ltd. All rights reserved.
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

(asdf:load-system :deoxybyte-systems)

(in-package :uk.co.deoxybyte-systems)

(asdf:defsystem genotype-call
  :name "genotype-call"
  :version "0.7.0"
  :author "Keith James"
  :licence "GPL v3"
  :depends-on (:deoxybyte-systems
               (:version :deoxybyte-io "0.12.0"))
  :in-order-to ((test-op (load-op :genotype-call :genotype-call-test)))
  :components ((:module :genotype-call
                        :serial t
                        :pathname "src/"
                        :components ((:file "package")
                                     (:file "utilities")
                                     (:file "gtc")
                                     (:file "bpm")
                                     (:file "sim")
                                     (:file "illuminus")
                                     (:file "intensity-processing")))
               (:lift-test-config :lift-tests
                                  :pathname "genotype-call-test"
                                  :target-system :genotype-call)))
