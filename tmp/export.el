(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(add-to-list 'package-archives '("org" . "https://orgmode.org/elpa/") t)
(package-initialize)
(when (not package-archive-contents)
  (package-refresh-contents))
;;
;; declare load-path
;; https://www.emacswiki.org/emacs/LoadPath
(let ((default-directory  "~/.emacs.d/lisp/"))
  (normal-top-level-add-to-load-path '("."))
  (normal-top-level-add-subdirs-to-load-path))
(let ((default-directory  "~/.emacs.d/elpa/"))
  (normal-top-level-add-to-load-path '("."))
  (normal-top-level-add-subdirs-to-load-path))
(let ((default-directory  "~/.emacs-packages/"))
  (normal-top-level-add-to-load-path '("."))
  (normal-top-level-add-subdirs-to-load-path))

;;

(require 'org)

;; Declare Babel languages
(org-babel-do-load-languages
 'org-babel-load-languages
 '(
   (ditaa . t)
   (dot .t)
   (emacs-lisp . t)
   (latex . t)
   (org . t)
   (python . t)
   (R . t)
   (shell . t)
   (sql .t)
   (sqlite . t)
   ))


(require 'snakemake-mode)

(defcustom snakemake-indent-field-offset nil
  "Offset for field indentation."
  :type 'integer)

(defcustom snakemake-indent-value-offset nil
  "Offset for field values that the line below the field key."
  :type 'integer)
