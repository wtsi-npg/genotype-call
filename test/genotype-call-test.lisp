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

(in-package :uk.ac.sanger.genotype-call-test)

(defparameter *example-snp-counts*
  '(("16" 4) ("2" 4) ("9" 2)))

(deftestsuite genotype-call-tests ()
  ())

(addtest (genotype-call-tests) read-bpm/1
  (with-open-file (stream (merge-pathnames "data/example.bpm.csv"))
    (let ((manifest (read-bpm stream))
          (chrs (mapcar #'first *example-snp-counts*))
          (counts (mapcar #'second *example-snp-counts*)))
      (ensure (= 10 (num-snps-of manifest)))
      (ensure (equal chrs (chromosomes-of manifest)))
      (ensure (equal counts (mapcar (lambda (chr)
                                      (num-snps-of manifest chr)) chrs))))))
