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

(defun generate-study (study-name study-size manifest-file
                       &key (x 1) (y 3) (genotype "AA") (basecall "AA")
                       (score 100.f0))
  "Writes GTC files for a fake study."
  (let ((bpm (load-bpm manifest-file :strict-ordering t))
        (meta-file (make-pathname :name study-name :type "json")))
    (with-open-file (meta meta-file :direction :output :if-exists :supersede
                          :if-does-not-exist :create)
      (loop
         for i from 0 below study-size
         for j = (num-snps-of bpm)
         for sample-name = (format nil "~a_~4,'0d" study-name i)
         for sample-urn = (format nil "urn:wtsi:~a" sample-name)
         for gtc-path = (make-pathname :name sample-name :type "gtc")
         collect gtc-path into gtc-paths
         collect (pairlis '(:result :uri :gender :gender-code)
                          (list (pathstring gtc-path) sample-urn
                                (if (oddp i)
                                    "Male"
                                    "Female")
                                (if (oddp i)
                                    1
                                    2)))
         into sample-meta
         do (write-gtc
             gtc-path (pathname-name manifest-file) sample-name
             (make-array (* 2 j) :element-type 'uint16 :initial-element x)
             (make-array (* 2 j) :element-type 'uint16 :initial-element y)
             (make-array j :initial-element genotype)
             (make-array j :initial-element basecall)
             (make-array j :initial-element score))
         finally (write-json-sample-specs sample-meta meta)
           (return (values gtc-paths meta-file))))))

(defun generate-manifest (filespec num-snps)
  "Writes a fake Beadpool Manifest for NUM-SNPS to FILESPEC."
  (flet ((make-name (i)
           (format nil "snp~7,'0d" i))
         (make-pos (chr i)
           (+ (* chr (expt 10 6)) i)))
    (with-open-file (bpm filespec :direction :output :if-exists :supersede
                         :if-does-not-exist :create)
      (write-line *bpm-header* bpm)
      (do* ((i 0 (1+ i))
            (per-chr (ceiling num-snps 23))
            (this-chr 0 (1+ this-chr))
            (snp-index 1 (1+ snp-index))
            (name (make-name snp-index) (make-name snp-index))
            (chromosome 1)
            (position (make-pos chromosome i) (make-pos chromosome i))
            (gentrain-score 1.f0)
            (snp "[G/A]")
            (ilmn-strand "TOP")
            (customer-strand "TOP")
            (norm-id 1))
           ((= num-snps i) filespec)
        (when (= per-chr this-chr)
          (setf this-chr 0
                chromosome (1+ chromosome)))
        (dolist (elt (intersperse (list snp-index
                                        name chromosome position
                                        gentrain-score snp
                                        ilmn-strand customer-strand norm-id)
                                  #\,))
          (princ elt bpm))
        (terpri bpm)))))

(defun write-gtc (filespec manifest-name sample-name x-intensities y-intensities
                  genotypes basecalls scores)
  "Writes a GTC file containing the supplied data."
  (check-arguments (stringp manifest-name) (manifest-name) "expected a string")
  (check-arguments (stringp sample-name) (sample-name) "expected a string")
  (check-arguments (every #'vectorp
                          (list x-intensities y-intensities genotypes
                                basecalls scores))
                   (x-intensities y-intensities genotypes basecalls scores) 
                   "expected all data in vectors")
  (let ((nx (length x-intensities))
        (ny (length y-intensities))
        (ng (length genotypes))
        (nb (length basecalls))
        (ns (length scores)))
    (check-arguments (= nx ny) (x-intensities y-intensities)
                     "intensity vectors must be the same length")
    (check-arguments (evenp nx) (x-intensities)
                     "there must be an even number of intensities")
    (let ((num-snps (/ nx 2)))
      (check-arguments (= ng num-snps) (genotypes)
                       "expected ~d genotypes, but found ~d" num-snps ng)
      (check-arguments (= nb num-snps) (basecalls)
                       "expected ~d basecalls, but found ~d" num-snps nb)
      (check-arguments (= ns num-snps) (scores)
                       "expected ~d scores, but found ~d" num-snps ns)
      (with-gtc (gtc filespec :direction :output :if-exists :supersede
                     :if-does-not-exist :create)
        (with-slots (stream)
            gtc
          (let* ((buffer (make-array 4 :element-type 'octet :initial-element 0))
                 (xform (make-xform 1 0.f0 0.f0 1.f0 1.f0)) ; no-op xform
                 (xforms (make-array 1 :initial-element xform))
                 (toc-entries '(:sample-name :snp-manifest
                                :normalization-xforms
                                :x-intensities :y-intensities
                                :genotypes :basecalls :genotype-scores))
                 (num-entries (1+ (length toc-entries))) ; +1 for :num-snps
                 (entry-size (+ 2 4))
                 (pad (+ 3 1 4 (* entry-size num-entries)))
                 (offsets (list 0
                                (1+ (length sample-name))   ; string
                                (1+ (length manifest-name)) ; string
                                (+ 4 (* 4 13))              ; single float
                                (+ 4 (* 2 nx))              ; uint16
                                (+ 4 (* 2 ny))              ; uint16
                                (+ 4 ng)                    ; byte
                                (+ 4 (* 2 nb)))))           ; byte pairs
            (write-magic *gtc-magic* stream)            ; 3 bytes of pad
            (write-version *gtc-version* stream buffer) ; 1 byte of pad
            (write-uint32 num-entries stream buffer)    ; 4 bytes of pad
            (write-gtc-toc-entry :num-snps num-snps stream buffer)
            (mapc (lambda (name offset)
                    (write-gtc-toc-entry name (+ pad offset) stream buffer))
                  toc-entries
                  (nreverse (maplist (lambda (x)
                                       (reduce #'+ x)) (nreverse offsets))))
            (write-gtc-string sample-name stream buffer)
            (write-gtc-string manifest-name stream buffer)
            (write-gtc-xforms xforms stream buffer)
            (write-gtc-intensities x-intensities stream buffer)
            (write-gtc-intensities y-intensities stream buffer)
            (write-gtc-genotypes genotypes stream buffer)
            (write-gtc-basecalls basecalls stream buffer)
            (write-gtc-genotype-scores scores stream buffer))))))
  filespec)
