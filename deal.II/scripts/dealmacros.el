;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; $Id$
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Emacs macros supporting deal.II programming
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun deal-indentation ()
"Install the indentation offsets used in the deal.II source code.
All contributions to the library should follow this indentation in
order to maintain a common look and feel as well as avoiding
unnecessary diffs in the archives.

The function is intended as a part of the cc-mode-startup-fun hook."

  (setq c-echo-semantic-information-p t)
  (setq c-basic-offset 2)
  (c-set-offset 'string                0)
  (c-set-offset 'defun-open            0)
  (c-set-offset 'defun-close           0)
  (c-set-offset 'defun-block-intro     '+)
  (c-set-offset 'class-open            0)
  (c-set-offset 'class-close           0)
  (c-set-offset 'inline-open           '+)
  (c-set-offset 'inline-close          0)
  (c-set-offset 'knr-argdecl-intro     '+)
  (c-set-offset 'knr-argdecl           0)
  (c-set-offset 'topmost-intro         0)
  (c-set-offset 'topmost-intro-cont    0)
  (c-set-offset 'member-init-intro     16)
  (c-set-offset 'member-init-cont      0)
  (c-set-offset 'inher-intro           '+)
  (c-set-offset 'inher-cont            'c-lineup-multi-inher)
  (c-set-offset 'block-close           0)
  (c-set-offset 'brace-list-open       0)
  (c-set-offset 'brace-list-close      0)
  (c-set-offset 'brace-list-intro      6)
  (c-set-offset 'brace-list-entry      0)
  (c-set-offset 'statement             'c-lineup-runin-statements)
  (c-set-offset 'statement-cont        'c-lineup-math)
  (c-set-offset 'statement-block-intro '+)
  (c-set-offset 'statement-case-intro  6)
  (c-set-offset 'statement-case-open   0)
  (c-set-offset 'substatement          '+)
  (c-set-offset 'substatement-open     '+)
  (c-set-offset 'case-label            '+)
  (c-set-offset 'access-label          '-)
  (c-set-offset 'label                 2)
  (c-set-offset 'do-while-closure      0)
  (c-set-offset 'else-clause           0)
  (c-set-offset 'arglist-intro         '+)
  (c-set-offset 'arglist-cont          0)
  (c-set-offset 'arglist-cont-nonempty 'c-lineup-arglist)
  (c-set-offset 'arglist-close         0)
  (c-set-offset 'stream-op             'c-lineup-streamop)
  (c-set-offset 'inclass               '++)
  (c-set-offset 'cpp-macro             -1000)
  (c-set-offset 'friend                0)

  (c-set-offset 'comment-intro         'c-lineup-comment)
  (c-set-offset 'c                     'c-lineup-C-comments)
  
  (c-set-offset 'objc-method-intro     -1000)
  (c-set-offset 'objc-method-args-cont 'c-lineup-ObjC-method-args)
  (c-set-offset 'objc-method-call-cont 'c-lineup-ObjC-method-call)
	
  (setq c-comment-only-line-offset 33)
  (setq c-hanging-comment-ender-p   nil)
  (setq c-hanging-comment-starter-p nil)

  (setq c-tab-always-indent t)
)

(defun deal-newline ()
"Setup emacs to automatically insert newlines and indentation before
and after semicolon or braces, like it is done in the deal.II coding style."
  (define-key c++-mode-map "\C-m" 'newline-and-indent)
  (setq c-hanging-braces-alist '((defun-open        . (before after))
                                 (defun-close       . (before after))
                                 (class-open        . (before after))
                                 (class-close       . (after))
                                 (inline-open       . (before after))
                                 (inline-close      . (before after))
                                 (block-open        . (before after))
                                 (block-close       . (after))
                                 (brace-list-open   . (nil))))
  (setq c-hanging-colons-alist '((member-init-intro . (after))
				 (inher-intro       . (after))
				 (label             . (after))
				 (case-label        . (after))
				 (access-label      . (after))))
  (setq c-cleanup-list '(empty-defun-braces
			 defun-close-semi
			 list-close-comma
			 scope-operator))
  (setq c-default-macroize-column 65)
)

(provide 'dealmacros)
