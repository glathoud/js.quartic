#!/usr/bin/env gsi-script

;;; Scheme - successfully tested with Gambit v4.2.8
;;; -*- mode:scheme; coding: utf-8 -*-
;;; Guillaume Lathoud - glathoud _at_ yahoo _dot_ fr
;;;
;;; Source code for ./index.html   
;;; The rendering of the mathematical formulas is done with MathJax,
;;; supposed to be present in a directory ../MathJax/
;;;
;;;
;;; Copyright 2013 Guillaume Lathoud
;;;
;;; Licensed under the Apache License, Version 2.0 (the "License");
;;; you may not use this file except in compliance with the License.
;;; You may obtain a copy of the License at
;;;
;;;     http://www.apache.org/licenses/LICENSE-2.0
;;;
;;; Unless required by applicable law or agreed to in writing, software
;;; distributed under the License is distributed on an "AS IS" BASIS,
;;; WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;;; See the License for the specific language governing permissions and
;;; limitations under the License.
;;;
;;; A copy of the Apache License Version 2.0 as of February 20th, 2013
;;; can be found in the file ./LICENSE.TXT

(include "../my_web.scm")

(define my_github "https://github.com/glathoud/js.quartic")

(define my_output_filename "index.html")

(define my_output_js_filename "index.js")
(with-output-to-file my_output_js_filename
  (lambda ()
    (display (file_concat_js 
              '("log.js"
                "complex.js"
                "solve_quartic.js"
                "solve_problem.js"
                "trackballcontrols.js"
                "three_view.js"
                )
              )
     )
    )
  )
(display my_output_js_filename)
(display "\n")

;(define my_title "Quartic equation resolution implementation, with a construction planning application")
(define my_title "Position the ramp of a construction site by solving a quartic equation")
(define my_keywords "math, quartic, construction")

(define my_license "
  This HTML page contains the mathematical description of a problem
  and its solution, an interactive numerical demo and 3D demo, and 
  links to the programs I wrote to implement the solution.


  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the \"License\");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an \"AS IS\" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
")

(define (my_pre_js . x)
  `( (pre (class "prettyprint lang-js")) ,@x)
)

(define (my_code_js . x)
  `( (code (class "prettyprint lang-js")) ,@x)
)

(define my_biblio_alist
   `(

     (quartic-history
      . (
         (ref . quartic-history)
         (title . ,(my_link 'quartic-history))
      )
     )

     (quartic-maximum
      . (
         (ref . quartic-maximum)
         (title . ,(my_link 'quartic-maximum))
      )
     )


     (quartic-wikipedia
      . (
         (ref . quartic-wikipedia)
         (title . ,(my_link 'quartic-wikipedia))
      )
     )

     (cubic-depressed-vieta
      . (
         (ref . cubic-depressed-vieta)
         (title . ,(my_link 'cubic-depressed-vieta))
      )
     )




     )
   )

(define my_local_bib (my_biblio my_biblio_alist))


(define def "\\buildrel\\triangle\\over =")

(define defspace (string-append "~" def "~"))


;;; Prepare the body

(define my_bodylist
  `( 

    
    ;; From https://github.com/blog/273-github-ribbons
    ( (a (href ,my_github)) 
      ( (img (style "position: fixed; top: 0; right: 0; border: 0; margin: 0; padding: 0;")
             (src "https://s3.amazonaws.com/github/ribbons/forkme_right_green_007200.png")
             (alt "Fork me on GitHub")
             )
        )
      )

    ;;

      ((div (class "topleftfixed topleftmenu"))

       (div ,(my_track "#main" "Top"))
       (div ,(my_track "#details" "Details"))
       (div ,(my_track "#horizproblem" "Geometry"))
       (div ,(my_track "#equation" "Equation"))
       (div ,(my_track "#solution" "Solution"))
       (div ,(my_track "#files" "Files"))
       (div "&nbsp;")
       (div ,(my_link_track 'home "&rArr; Back home"))
       (div ,(my_link_track 'email "&rArr; email"))
;       (div "&nbsp;")
;       (div ,my_flatter)
       )
 
      ( (div (id main))

      (h2 ,my_title)
      ( (h3 (class fr)) "by " ,my_author ", February 2013")
      ( (div (class clear)))

      ; Tried ( (a (href "javascript:(function(){var script=document.createElement('script');script.type='text/javascript';script.src='https://github.com/zz85/zz85-bookmarklets/raw/master/js/ThreeInspector.js';document.body.appendChild(script);})()")) "Three.js Scene Inspector")


      ( (div (id three-view)))

      ( (h2 (id problem)) "Problem")

      ( (form (name problem-input))

      (p "Determine the position of the blue ramp with:"

         (blockquote "width&nbsp;" ( (input (type number) (class number-input) (name blue-width) (value 3.5) (min 0) (max 1e4) (step 0.1))) "&nbsp;meters and slope&nbsp;" ( (input (type number) (class number-input) (name blue-slope) (value 25) (min 1) (max 1e5))) "&nbsp;%")
         
         "when the red slopes have:"
         
         (blockquote "depth " ( (input (type number) (class number-input) (name red-depth) (value 3.5) (min 0) (max 1e4) (step 0.1))) "&nbsp;meters and slope " ( (input (type number) (class number-input) (name red-slope) (value 60) (min 1) (max 89))) "&nbsp;degrees.")
         
         )
      )

      ( (h2 (id result)) "Result")
      
      ( (form (name result))
        
        (p 
         "The blue ramp has:"
         (blockquote "length " ( (input (type text) (readonly readonly) (class number-output) (name blue-length)))
                     " meters and angle " ( (input (type text) (readonly readonly) (class number-output) (name blue-angle))) 
                     " degrees."
                     )
         " and the small side of the green triangle has:"
         (blockquote "length " ( (input (type text) (readonly readonly) (class number-output) (name green-small-side-length)))
                     " meters."
                     )
         )
        )

      ( (h2 (id details)) "Details" )

      (p "Below, the " ( (a (href "#problem")) "problem") " is described geometrically in the
" ( (a (href "#horizproblem") ) "horizontal plane" ) ", turned into a 4th power polynomial " ( (a (href "#equation")) "equation" ) ", of which a
direct, " ( (a (href "#solution")) "general solution") " is known since the
Renaissance works of Ferrari (and Vieta).")
      
      (p ( (a (href "#files")) "Files") " are then shortly described, which implement this general solution.")

      ( (h3 (id problem-details)) "Problem")

      (p "The goal is to determine the position of blue ramp going down from the top of the red slopes to their bottom (all three are rectangles). ")


      (p "Both steep red slopes have the same  \\( depth\\) and \\(slope_{red}\\), thus the same projected «horizontal width»: \\[
 qq_{red}" ,defspace " \\frac{depth}{slope_{red}}
\\] where for example \\(slope_{red} = tan(\\pi\\frac{60}{180})\\) (60-degree red slopes).")

      (p "The blue ramp has a given \\(width_{blue}\\) and a given \\( slope_{blue}\\), and thus a projected «horizontal length»: \\[ 
horizlength_{blue}" ,defspace "\\frac{depth}{slope_{blue}}
\\] where \\(slope_{blue} < slope_{red} \\), for example \\(slope_{blue} = 0.25 ~~(25\\%)\\).")
      
      ( (h3 (id horizproblem)) "Problem, projected onto the horizontal plane")
      
      (div ( (img (src horiz.jpg))))

      (h4 "Constraints")

      (p "The green triangle is horizontal. The leftmost corner of the blue ramp \\(i \\cdot zz\\) must be at the top of the red slope.")

      (p "The origin \\(O\\), the nearby corner of the blue ramp \\( yy \\cdot exp(i \\cdot angle) \\), and the rightmost corner of the blue ramp \\( xx + i \\cdot qq \\) must be colinear. The latter must also be at the bottom of the red slope.")

;      (p "Since the ramp is rectangular, the remaining corner is: \\( i \\cdot zz + horizlength \\cdot exp(i \\cdot angle) \\).")

      (p "The problem is fully determined 
by \\(qq\\), \\(horizlength\\) and \\(width\\). To&nbsp;find the position of the blue ramp, we need to determine \\(xx\\), \\(yy\\), \\(zz\\) and \\(angle\\).")

            
      ( (h3 (id equation)) "Equation")

      (p "Since the ramp is rectangular, \\( angle \\) appears at several locations. Expressing for example \\(tan(angle)\\), one can show that: \\[ qq \\cdot width=xx \\cdot yy \\]
      Squaring both sides and using the Pythagorean theorem to express \\(xx^2\\), one obtains:
 
      \\[ f(yy)=0 \\]

      where: \\[ f(yy)" ,defspace "yy^2 \\cdot ((horizlength+yy)^2-qq^2)-qq^2 \\cdot width^2 \\]

which is a 4th power polynomial in \\( yy \\). The other variables \\(xx\\), \\(zz\\) and \\(angle\\) 
can all be derived from \\( yy \\).
")
      
      ( (h3 (id solution)) "Direct solution")

      (p "Fortunately, the 4th power is the highest degree that can be 
algebraically solved&nbsp;"   ,(my_local_bib 'quartic-history) ,(my_local_bib 'quartic-maximum) " using a direct, general solution, which I decided 
to implement, so as to have reusable code. I chose Ferrari's solution&nbsp;" ,(my_local_bib 'quartic-wikipedia) ", including Vieta's substitution in the depressed cubic case&nbsp;" ,(my_local_bib 'cubic-depressed-vieta) ".")

      ( (h3 (id files)) "Files")

      (p "Code I wrote to solve the problem:")
      (ul
       (li ( (a (href "log.js")) "log.js") " and " ( (a (href "complex.js")) "complex.js") ": Base code.")
       (li ( (a (href "solve_quartic.js")) "solve_quartic.js" ) ": " (strong "General solution of the quartic equation."))
       (li ( (a (href "solve_problem.js")) "solve_problem.js" ) ": Application to the specific ramp problem.")
       (li ( (a (href "./LICENSE.TXT")) "LICENSE.TXT") ": Apache 2.0 License." )
       )

      (p "Third-party software, used to render this article:")
      (ul
       (li ( (a (href "three.js")) "three.js" ) " and " ( (a (href "three.min.js")) "three.min.js") ": 3D rendering engine.")
       (li ( (a (href "trackballcontrols.js")) "trackballcontrols.js") ": Extension to rotate the 3D model.")
       (li ( (a (href "../MathJax/MathJax.js")) "MathJax.js") ": Mathematical formula rendering engine.")
       )


      ( (h3 (id ack)) "Acknowledgments")

      (p "Thanks to Felix Schädler for the initial request, and to the authors of " ,(my_link_track 'three_js "three.js") " and " ,(my_link_track 'mathjax "MathJax") " for their excellent softwares." )

      ( (h3 (id ref)) "References")

     ,@(my_local_bib)

     ( (h3 (id license)) "License" )

     (pre ,my_license )

     (p "File " ( (a (href "LICENSE.TXT")) "./LICENSE.TXT"))
)

))



;;; Javascript part

(define my_head_tail_list
  `(
    ( (script (type text/javascript) (src "../prettify/prettify.js") ) )
    ( (link (href "../prettify/prettify.css") (type text/css) (rel stylesheet)) )
    "
<!--

"
,my_license 
"
-->
"
    )
  )

(define my_body_tail_list
  `( 
    ( (script (type text/javascript)) "setTimeout( prettyPrint );" ) 
    ( (script (type text/javascript) (src three.min.js)))
    ( (script (type text/javascript) (src ../MathJax/MathJax.js?config=default)))
    ( (script (type text/javascript) (src index.js)))
    )
)

;;; Prepare the contents: an alist

(define my_contents
  `( 
    (is_html   . #t)
    (title     . ,my_title)
    (author    . ,my_author)
    (email     . ,my_email)
    (keywords  . ,my_keywords)
    (head_css  . "body { margin: auto; width: 600px; }
@media only screen and (max-width: 1040px) { body { margin-left: 220px; } }
#main { width: 100%; margin: auto; }
@media only screen { .bib-view-single { display: none; } }

#three-view { width: 100%; height: 400px; cursor: pointer; }

.number-input, .number-output { width: 45px; }
pre, body { color: #232; }
")
    (head_tail_list . ,my_head_tail_list)
    (bodyprop  . ((class "")))   ; className "loading" if needed
    (bodylist  . ,my_bodylist)
    (body_tail_list . ,my_body_tail_list)
    )
  )

;;; Feed the contents into the layout and display the resulting XHTML string

(define result (my_layout my_contents))

(with-output-to-file my_output_filename
  (lambda ()
    (display result)
    )
  )

(display (string-append "Wrote: " my_output_filename "\n"))
