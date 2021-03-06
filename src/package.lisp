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

(in-package :cl-user)

(defpackage :uk.ac.sanger.genotype-call
  (:nicknames #:genotype-call)
  (:use #:common-lisp #:deoxybyte-utilities #:deoxybyte-io)
  (:export
   
   #:*xform-fields*

   ;; CLI
   #:cli
   
   ;; BPM
   #:bpm
   #:load-bpm
   #:read-bpm
   #:alleles-of
   #:chromosomes-of
   #:snps-of
   #:num-snps-of
   #:valid-alleles-p
   #:has-chromosome-p
   #:chromosome-boundaries
   #:location<

   #:make-chromosome-p

   ;; SNP
   #:make-snp
   #:snp
   #:snp-index
   #:snp-name
   #:snp-chromosome
   #:snp-position
   #:snp-alleles
   #:snp-ilmn-strand
   #:snp-cust-strand
   #:snp-norm-id
   #:snp-norm-rank

   ;; GTC
   #:gtc
   #:with-gtc
   #:gtc-open
   #:gtc-close
   #:read-gtc
   #:data-field-of
   #:toc-of
   #:normalize
   #:xform
   #:xform-p

   ;; SIM
   #:sim
   #:with-sim
   #:sim-open
   #:sim-close
   #:name-size-of
   #:num-samples-of
   #:num-probes-of
   #:num-channels-of
   #:format-of
   #:read-intensities

   ;; BED
   #:bed
   #:with-bed
   #:bed-open
   #:bed-close
   #:orientation-of

   ;; Intensity processing
   #:copy-intensities
   #:gtc-to-sim
   #:sim-to-illuminus
   #:sim-to-genosnp
   ;; Genotype processing
   #:copy-genotypes
   #:gtc-to-bed
   #:plink-pathname

   ;; Utilities
   #:version-of)
  (:documentation "This package provides functions for reading
  Genotype Call (GTC) files created by Illumina's AutoCall
  software."))
