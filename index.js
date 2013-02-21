//#---------- file: "log.js"

/*
  Function that prints log messages.
  

  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

(function (G) {

    if (G.log)
        return;

    G.log = log;

    function log(s)
    {
        if (typeof console !== 'undefined')
            console.log(s);  // Browsers
        else
            print(s);  // V8, Rhino
    }

})(this);

//#---------- end of file: "log.js"

//#---------- file: "complex.js"

/*
  Functions to create and manipulate complex numbers.
  

  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

(function (/*global*/G) {

    if (G.Complex)
        return;

    var epsilon = 1e-5;

    // Public API

    G.cplx = cplx;
    G.Complex = Complex;

    G.cadd  = cadd;
    G.csub  = csub;
    G.cmul  = cmul;
    G.cdiv  = cdiv;
    G.cpol  = cpol;
    G.croot = croot;

    // Implementation

    function cadd(a,b) { return cplx(a).add( cplx(b) ); } 
    function csub(a,b) { return cplx(a).sub( cplx(b) ); } 
    function cmul(a,b) { return cplx(a).mul( cplx(b) ); } 
    function cdiv(a,b) { return cplx(a).div( cplx(b) ); } 

    function croot(a, n, /*?optional 0<k<=n?*/k) { return cplx(a).root(n, k); } 

    function cplx(a,b) { return a instanceof Complex ? a : new Complex(a, b); } 

    function cpol(r,angle) { return new Complex(r*Math.cos(angle),r*Math.sin(angle)); } 

    function Complex(re,im) 
    {
        this.re = re || 0;
        this.im = im || 0;
    }

    var Cp = Complex.prototype;
    Cp.toString = Cp.toValue = function () { 
        throw new Error('Complex does not support toString / toValue - to ease debugging'); 
    };
    Cp.str = function () { 
        var re = this.re
        ,   im = this.im
        ;
        return this.ltEps()  ?  '<ltEps: ' + re.toPrecision(4) + '+i*' + im.toPrecision(4) + '>'  :  re + '+i*' + im; 
    };
    Cp.add = function (other) {
        return cplx(this.re+other.re,this.im+other.im); 
    };
    Cp.sub = function (other) { 
        return cplx(this.re-other.re,this.im-other.im); 
    };
    Cp.mul = function (other) { 
        return cplx(this.re*other.re-this.im*other.im, 
                    this.re*other.im+this.im*other.re
                   ); 
    };
    Cp.div = function (other) { 
        return this.mul(other.conj()).mul( cplx( 1 / other.r2() ) );
    };
    Cp.conj = function () {
        return cplx(this.re, -this.im);
    }
    Cp.root = function (/*integer*/n, /*?integer?*/k) {
        // Just return *one* of the roots
        k != null  ||  (k = 0); // Default value

        if (typeof n !== 'number'  ||  n !== ~~n)
            throw new Error('complex root: n must be an integer. You gave: ' + n);
        if (typeof k !== 'number'  ||  k !== ~~k)
            throw new Error('complex root: k must be an integer. You gave: ' + k);


        
        var   r = Math.pow(this.r2(), 1/(2*n))
        , angle = Math.atan2(this.im,this.re) / n
        ;
        return cpol(r,angle + k * 2 * Math.PI / n );
    }
    Cp.r2 = function ()  {
        return this.re*this.re+this.im*this.im;
    }
    Cp.ltEps = function () {
        return this.r2() < epsilon*epsilon;
    }
    Cp.isReal = function () {
        return this.ltEps()  ||  ( this.im * this.im / this.r2() ) < epsilon * epsilon;
    }

})(this);

//#---------- end of file: "complex.js"

//#---------- file: "solve_quartic.js"

/*
  ECMAScript implementation of the algebraic resolution of quartic
  equations, based on the Renaissance works of Ferrari and Vieta.

  
  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

/*global load Quartic cadd cmul cdiv csub croot log cplx*/
if (typeof load !== 'undefined')
{
    load('log.js');
    load('complex.js');
}

var quartic;
(function (G) {

    if (quartic)
        return;

    // API

    quartic = {
        getRootsQuartic          : getRootsQuartic
        , getRootsCubic          : getRootsCubic
        , getRootsDepressedCubic : getRootsDepressedCubic
        , getRootsQuadratic      : getRootsQuadratic
    };
    
    // API Implementation

    function getRootsQuartic( /*number | Complex*/A, B, C, D, E, /*?boolean?*/verbose )
    // Find the roots of a quartic polynomial 
    //
    // f(x) = A x^4 + B x^3 + C x^2 + D x + E
    //
    // that is the four values x1, x2, x3, x4
    // verifying f(x) = 0
    // 
    // The implementation
    // below is general (complexes all the way).
    //
    // It uses Ferrari's approach:
    // http://en.wikipedia.org/wiki/Quartic_function
    // 
    // including Vieta's substitution for the depressed cubic:
    // http://en.wikipedia.org/wiki/Cubic_equation#Vieta.27s_substitution
    // 
    // Guillaume Lathoud
    // 2012-12-17
    {
        // First, convert to a depressed quartic:
        //
        //   x = u - B / 4A
        // 
        //   u^4 + alpha * u^2 + beta * u + gamma = 0
        
        var A2 = cmul( A, A )
        , A3 = cmul( A, A2 )
        , A4 = cmul( A, A3 )

        , B2 = cmul( B, B )
        , B3 = cmul( B, B2 )
        , B4 = cmul( B, B3 )

        , alpha = cmul( -3, B2 ).div( cmul( 8, A2) ).add( C.div( A ) )
        , beta  = B3.div( cmul( 8, A3 ) ).sub( B.mul(C).div( cmul( 2, A2 ))).add(D.div(A))
        , gamma = cmul( -3, B4 ).div( cmul( 256, A4 ) ).add( C.mul(B2).div( cmul(16,A3))).sub(B.mul(D).div(cmul(4,A2))).add(E.div(A))
        ;
        
        if (verbose)
        {
            log('');
            log('depressed quartic:');
            log('alpha ' + alpha.str());
            log('beta ' + beta.str());
            log('gamma ' + gamma.str());
        }

        var u1,u2,u3,u4;  // Where the `u` roots values will be stored

        if (beta.ltEps())
            quartic_solve_biquadratic();
        else if (gamma.ltEps())
            quartic_solve_depressed_cubic();
        else
            quartic_solve_non_degenerate_quartic();

        var bdfa = B.div( cmul( 4, A ));
        if (verbose)
        {
            log('');
            log('B/4A: ' + bdfa.str());
        }          
        
        var x1 = u1.sub(bdfa)
        ,   x2 = u2.sub(bdfa)
        ,   x3 = u3.sub(bdfa)
        ,   x4 = u4.sub(bdfa)
        ;
        
        return [ x1, x2, x3, x4 ];

        function quartic_solve_depressed_cubic()
        // u4 = 0 is a root.
        //
        // For the 3 other roots:
        // 
        //   u^3 + alpha*u + beta = 0      the depressed cubic equation
        {
            if (verbose)
            {
                log('');
                log('quartic_solve_depressed_cubic');
            }
            
            var uroots = getRootsDepressedCubic( 1, alpha, beta, verbose );
            u1 = uroots[ 0 ];
            u2 = uroots[ 1 ];
            u3 = uroots[ 2 ];
            u4 = cplx( 0 );
        }

        function quartic_solve_biquadratic()
        // (u^2)^2 + alpha * (u^2) + gamma = 0
        {
            if (verbose)
            {
                log('');
                log('quartic_solve_biquadratic');
            }
            var qroots = getRootsQuadratic( 1, alpha, gamma, verbose )
            ,     u2_1 = qroots[ 0 ]
            ,     u2_2 = qroots[ 1 ]
            ;
            u1 = biquadratic_solution( u2_1, -1 );
            u2 = biquadratic_solution( u2_1, +1 );
            u3 = biquadratic_solution( u2_2, -1 );
            u4 = biquadratic_solution( u2_2, +1 );
        }

        function biquadratic_solution( u2, sign )
        {
            return cmul( sign > 0 ? +1: -1 , croot( u2, 2 ) );
        }

        function quartic_solve_non_degenerate_quartic()
        {
            // nested depressed cubic:

            var alpha2 = alpha.mul(alpha)
            ,   alpha3 = alpha.mul(alpha2)
            
            ,   P = cdiv( alpha2, -12 ).sub( gamma )
            ,   Q = cdiv( alpha3, -108 ).add( cdiv( alpha.mul(gamma), 3 )).sub( cdiv( beta.mul(beta), 8 ) )
            
            ,   tmp = croot( cdiv( Q.mul(Q), 4 ).add( cdiv( P.mul(P).mul(P), 27 ))
                             , 2 
                           )
            ,   R1 = cadd( cdiv( Q, -2 ), tmp)
            //  R2 = csub( cdiv( Q, -2 ), tmp);  // R1 enough

            ,   U = croot(R1, 3)
            
            ,   V = U.ltEps() ? croot(Q,3) : cdiv( cdiv(P,3),U)
            
            ,   y = cmul( -5/6, alpha ).add(U).sub( V )
            ;
            if (verbose)
            {
                log('');
                log('nested depressed cubic:');
                log('P ' + P.str());
                log('Q ' + Q.str());
                log('U ' + U.str());
                log('V ' + V.str());
                log('y ' + y.str());
            }

            // Folding the second perfect square

            var    W = croot(cadd(alpha,cmul(2,y)), 2)
            ,   taty = cadd( cmul(3,alpha) , cmul(2, y) )
            ,   tbw  = cdiv( cmul(2,beta) , W )
            ;
            if (verbose)
            {
                log('');
                log('Folding the second perfect square');
                log('W ' + W.str());
            }
            u1 = Wbtt_2_u( W, bdfa, taty, tbw, -1, -1 );
            u2 = Wbtt_2_u( W, bdfa, taty, tbw, -1, +1 );
            u3 = Wbtt_2_u( W, bdfa, taty, tbw, +1, -1 );
            u4 = Wbtt_2_u( W, bdfa, taty, tbw, +1, +1 );

            function Wbtt_2_u( W, bdfa, taty, tbw, s, t )
            {
                return cmul
                ( 0.5
                  , cadd
                  ( cmul( s, W )
                    , cmul
                    ( t
                      , croot
                      ( cmul
                        ( -1
                          , cadd( taty, cmul( s, tbw ) )
                        )
                        , 2
                      )
                    )
                  )
                );
            }
        }

    } // end of: getRootsQuartic

    function getRootsCubic( A, B, C, D, verbose )
    // Find the three roots of the cubic equation:
    // 
    //   A * x^3 + B * x^2 + C * x + D = 0
    //
    // Let:  x = u - B / 3A
    //
    // Then:
    // 
    //   A * u^3 + P * u + Q = 0
    // 
    // http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    {
        var A2 = cmul( A, A )
        ,   B2 = cmul( B, B )
        ,   B3 = B2.mul( B )
        
        ,   P  = csub( C, B2.div( cmul( 3, A ) ) )
        
        ,   Q  = cmul( 2, B3 )
            .sub( cmul( 9, cmul( A, cmul( B, C ) ) ) )
            .add( cmul( 27, cmul( A2, D ) ) )
            .div( cmul( 27, A2 ) )
        
        , dcroots = getRootsDepressedCubic( A, P, Q, verbose )
        ,      u1 = dcroots[ 0 ]
        ,      u2 = dcroots[ 1 ]
        ,      u3 = dcroots[ 2 ]

        ,    bdta = cdiv( B, cmul( 3, A ) )
        ,      x1 = u1.sub( bdta )
        ,      x2 = u2.sub( bdta )
        ,      x3 = u3.sub( bdta )
        ;  
        if (verbose)
        {
            log('');
            log('bdta: ' + bdta.str());
        }
        return [ x1, x2, x3 ];
    }


    function getRootsDepressedCubic( A, C, D, verbose )
    // Find the three roots of the depressed cubic equation:
    // 
    //   A * x^3 + C * x + D = 0
    // 
    // Let:
    //   alpha = C / A  
    //   beta  = D / A
    //
    // then:
    //   x^3 + alpha * x + beta = 0
    // 
    // The Vieta's substitution [vieta]:
    // 
    //   x = w - alpha / 3w
    //
    // results in the equation:
    //
    //   w^3 + beta - alpha^3 / (27 * w^3) = 0
    // 
    //   (w^3)^2 + beta (w^3) - alpha^3 / 27 = 0
    // 
    // -> solve the quadratic to find one value of w^3
    // -> compute the 3 roots w1 w2 w3 and thus x1 x2 x3
    //
    // [vieta] http://en.wikipedia.org/wiki/Cubic_equation#Vieta.27s_substitution
    {
        var alpha = cdiv( C, A )
        ,   beta  = cdiv( D, A )
        ,   delta = quadratic_delta( 1
                                     , beta
                                     , cdiv( alpha.mul(alpha).mul(alpha), -27)
                                   )
        , wcube = quadratic_root(1, beta, delta, +1)
        , w1 = wcube.root(3,0)
        , w2 = wcube.root(3,1)
        , w3 = wcube.root(3,2)
        , x1 = w2u( w1 )
        , x2 = w2u( w2 )
        , x3 = w2u( w3 )
        ;

        if (verbose)
        {
            log('');
            log('wcube ' + wcube.str());
            log('w1 ' + w1.str() + ' ' + (w1.mul(w1).mul(w1).sub(wcube)).ltEps());
            log('w2 ' + w2.str() + ' ' + (w2.mul(w2).mul(w2).sub(wcube)).ltEps())
            log('w3 ' + w3.str() + ' ' + (w3.mul(w3).mul(w3).sub(wcube)).ltEps());
        }
        
        return [x1, x2, x3];
        
        function w2u( w )
        {
            return w.sub( alpha.div( cmul(3,w) ) );
        }
    }


    function getRootsQuadratic( A, B, C, verbose )
    // A * x^2 + B * x + C = 0
    {
        var delta = quadratic_delta( A, B, C )
        ,   x1    = quadratic_root( A, B, delta, +1 )
        ,   x2    = quadratic_root( A, B, delta, -1 )
        ;
        if (verbose)
        {
            log('');
            log('delta: ' + delta.str() );
            log('x1:    ' + x1.str() );
            log('x2:    ' + x2.str() );
        }
        return [ x1, x2 ];
    }
    
    // ----------------------------------------------------------------------
    // Private details

    function quadratic_delta(a,b,c)
    {
        return cmul( b, b ).sub( cmul( 4, cmul( a, c ) ) )
    }

    function quadratic_root( a,b,delta,sign )
    // a*x^2 + b*x + c = 0
    // delta = b^2-4*a*c
    {
        return cdiv(
            cmul(-1, b)[ sign > 0 ? 'add' : 'sub']( croot( delta, 2 ))
            , cmul(2, a)
        );
    }

})();


//#---------- end of file: "solve_quartic.js"

//#---------- file: "solve_problem.js"

/*
  Solution to the specific ramp problem described in ./index.html
  Required files:  ./quartic.js  ./complex.js  ./log.js
  

  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

/*global quartic try_to_solve*/

function try_to_solve( /*object*/problem )
{
    var horizlength = problem.horizlength
    ,   width  = problem.width
    ,   qq     = problem.qq
    ;
    // Make sure we have all numbers
    horizlength.toPrecision.call.a , width.toPrecision.call.a , qq.toPrecision.call.a;

    /*

      Problem, viewed from above, projected onto the horizontal plane:

      |<---qq-->|
      | steep slope
      | (60deg) |               i*zz + horizlength * exp(i*angle)
      |         |
      i*zz      |          rampe 
      |         |  rampe
      |     rampe_(slope not steep)_(xx+i*qq)__________
      |                                                ^
      |   yy*exp(i*angle)         steep slope (60deg)  qq
      o________________________________________________v 

      Goal: Build a *rectangular* rampe of not-steep-slope 25% to go down
      some depth.  The above graph represents the "steep slope" and the four
      points of the rampe, projected in the horizontal plane.

      Projected onto the horizontal complex plane, the rampe is a rectangle
      `horizlength*with` with the four edge points:

      i*zz
      yy*exp(i*angle)
      xx+i*qq
      i*zz + horizlength * exp(i*angle)

      Assumptions/goals: the three points o (origin), yy*exp(i*angle) and
      xx+i*qq must be aligned. Also, i*zz coincides with the top of the
      steep slope and xx+i*qq coincides with the bottom of the steep slope.

      Note that horizlength must be > qq so that the
      slope of the rampe is less than the slope of 60deg.

      One can show that qq*width=xx*yy  (look e.g. at tan(angle))
      and that:

      f(yy)=0

      where f(yy):=yy^2*((horizlength+yy)^2-qq^2)-qq*qq*width*width

      The problem was given to me initially by Felix Schaedler.
      I solve f(yy)=0 first quickly, using dichotomy (see ./OLD/).

      Below, a direct solution that implements Ferrari's approach to find
      the 4 roots of any quartic equation (complex or real coefficients - in
      our application we have real coefficients).

      Ferrari's approach: http://en.wikipedia.org/wiki/Quartic_function
      The implementation below is general (complexes all the way).

      Guillaume Lathoud
      2012-11-28

      Numerical results:

      application-specific: look for a positive real root
      yy_solution 0.49281701552540547
      f(yy_solution) 0+i*0
      f(yy_solution).ltEps() true

      angle_solution_radians 0.13988524363241306
      angle_solution_degrees 8.014834076296541
      xx_solution 14.35125122462075
      xx_solution-qq 12.33052528245706
      zz_solution 3.5345252313134456
    */

    // Use case values A B C D E are fed into the generic complex
    // implementation, which finds the four roots x1 x2 x3 x4 of the
    // quartic equation f(yy)=0
    
    
    // The problem described above (rampe)
    // Application-specific values (reals)      
    //
    // This leads to a non-degenerate quartic.
    
    A = cplx( 1 );
    B = cplx( 2 * horizlength );
    C = cplx( horizlength*horizlength - qq*qq );
    D = cplx( 0 );
    E = cplx( -qq*qq*width*width );
    
    function f(yy)
    {
        yy = cplx(yy);
        var yy2 = yy.mul(yy)
        ,   yy3 = yy.mul(yy2)
        ,   yy4 = yy.mul(yy3)
        ;
        return cmul( A, yy4 ).add( cmul( B, yy3 ) ).add( cmul( C, yy2 ) ).add( cmul( D, yy ) ).add( E );
    }

    // Solve: find four roots

    var roots = quartic.getRootsQuartic( A, B, C, D, E, /*verbose:*/true )
    ,   x1    = roots[ 0 ]
    ,   x2    = roots[ 1 ]
    ,   x3    = roots[ 2 ]
    ,   x4    = roots[ 3 ]
    ;

    // Check that these really are roots of f(yy)

    var fx1 = f(x1), fx2 = f(x2), fx3 = f(x3), fx4 = f(x4)
    , verif1 = fx1.ltEps()
    , verif2 = fx2.ltEps()
    , verif3 = fx3.ltEps()
    , verif4 = fx4.ltEps()
    ;

    log('');
    log('roots:');
    log('x1 ' + x1.str() + ' -> f(x1) ' + fx1.str() + ' verif1:' + verif1);
    log('x2 ' + x2.str() + ' -> f(x2) ' + fx2.str() + ' verif2:' + verif2);
    log('x3 ' + x3.str() + ' -> f(x3) ' + fx3.str() + ' verif3:' + verif3);
    log('x4 ' + x4.str() + ' -> f(x4) ' + fx4.str() + ' verif4:' + verif4);

    if (!(verif1*verif2*verif3*verif4))
        throw new Error('At least one of x1,x2,x3,x4 is not a root!');

    log('');
    log('--> ROOTS OK!');
    log('');


    // Application-specific: try to look for a positive real root

    var yy_solution = null
    , candidates = [ x1, x2, x3, x4 ]
    ;
    for (var i = candidates.length; i--;)
    {
        var c = candidates[i];
        if (!c.isReal())
            continue;

        yy_solution = Math.max( yy_solution || -Infinity, c.re );
    }

    log('');
    log('application-specific: try to look for a positive real root');
    log('yy_solution ' + yy_solution);

    var solution = null;
    if (yy_solution == null)
        log('-> not found.');
    else
    {
        log('-> found.');
        log('');
        log('f(yy_solution) ' + f(yy_solution).str());
        
        var verification = f(yy_solution).ltEps();
        log('f(yy_solution).ltEps() ' + verification);
        if(!verification)
            throw new Error('Something is wrong!');

        log('');

        // Application-specific: derive other values

        angle_solution_radians = Math.atan2(yy_solution,width);
        angle_solution_degrees = angle_solution_radians / Math.PI * 180;
        xx_solution = qq*width/yy_solution;
        zz_solution = Math.sqrt(width*width+yy_solution*yy_solution);

        log('angle_solution_radians ' + angle_solution_radians);
        log('angle_solution_degrees ' + angle_solution_degrees);
        log('xx_solution ' + xx_solution);
        log('xx_solution-qq ' + (xx_solution-qq));
        log('zz_solution ' + zz_solution);

        solution = {
            xx   : xx_solution
            , yy : yy_solution
            , zz : zz_solution
            , angle_radians : angle_solution_radians
            , angle_degrees : angle_solution_degrees
        };
        
    }

    log('');  

    return solution;
}

//#---------- end of file: "solve_problem.js"

//#---------- file: "trackballcontrols.js"

/**
 * @author Eberhard Graether / http://egraether.com/
 */

THREE.TrackballControls = function ( object, domElement ) {

	THREE.EventDispatcher.call( this );

	var _this = this;
	var STATE = { NONE: -1, ROTATE: 0, ZOOM: 1, PAN: 2, TOUCH_ROTATE: 3, TOUCH_ZOOM: 4, TOUCH_PAN: 5 };

	this.object = object;
	this.domElement = ( domElement !== undefined ) ? domElement : document;

	// API

	this.enabled = true;

	this.screen = { width: 0, height: 0, offsetLeft: 0, offsetTop: 0 };
	this.radius = ( this.screen.width + this.screen.height ) / 4;

	this.rotateSpeed = 1.0;
	this.zoomSpeed = 1.2;
	this.panSpeed = 0.3;

	this.noRotate = false;
	this.noZoom = false;
	this.noPan = false;

	this.staticMoving = false;
	this.dynamicDampingFactor = 0.2;

	this.minDistance = 0;
	this.maxDistance = Infinity;

	this.keys = [ 65 /*A*/, 83 /*S*/, 68 /*D*/ ];

	// internals

	this.target = new THREE.Vector3();

	var lastPosition = new THREE.Vector3();

	var _state = STATE.NONE,
	_prevState = STATE.NONE,

	_eye = new THREE.Vector3(),

	_rotateStart = new THREE.Vector3(),
	_rotateEnd = new THREE.Vector3(),

	_zoomStart = new THREE.Vector2(),
	_zoomEnd = new THREE.Vector2(),

	_touchZoomDistanceStart = 0,
	_touchZoomDistanceEnd = 0,

	_panStart = new THREE.Vector2(),
	_panEnd = new THREE.Vector2();

	// events

	var changeEvent = { type: 'change' };


	// methods

	this.handleResize = function () {

		this.screen.width = window.innerWidth;
		this.screen.height = window.innerHeight;

		this.screen.offsetLeft = 0;
		this.screen.offsetTop = 0;

		this.radius = ( this.screen.width + this.screen.height ) / 4;

	};

	this.handleEvent = function ( event ) {

		if ( typeof this[ event.type ] == 'function' ) {

			this[ event.type ]( event );

		}

	};

	this.getMouseOnScreen = function ( clientX, clientY ) {

		return new THREE.Vector2(
			( clientX - _this.screen.offsetLeft ) / _this.radius * 0.5,
			( clientY - _this.screen.offsetTop ) / _this.radius * 0.5
		);

	};

	this.getMouseProjectionOnBall = function ( clientX, clientY ) {

		var mouseOnBall = new THREE.Vector3(
			( clientX - _this.screen.width * 0.5 - _this.screen.offsetLeft ) / _this.radius,
			( _this.screen.height * 0.5 + _this.screen.offsetTop - clientY ) / _this.radius,
			0.0
		);

		var length = mouseOnBall.length();

		if ( length > 1.0 ) {

			mouseOnBall.normalize();

		} else {

			mouseOnBall.z = Math.sqrt( 1.0 - length * length );

		}

		_eye.copy( _this.object.position ).sub( _this.target );

		var projection = _this.object.up.clone().setLength( mouseOnBall.y );
		projection.add( _this.object.up.clone().cross( _eye ).setLength( mouseOnBall.x ) );
		projection.add( _eye.setLength( mouseOnBall.z ) );

		return projection;

	};

	this.rotateCamera = function () {

		var angle = Math.acos( _rotateStart.dot( _rotateEnd ) / _rotateStart.length() / _rotateEnd.length() );

		if ( angle ) {

			var axis = ( new THREE.Vector3() ).crossVectors( _rotateStart, _rotateEnd ).normalize(),
				quaternion = new THREE.Quaternion();

			angle *= _this.rotateSpeed;

			quaternion.setFromAxisAngle( axis, -angle );

			_eye.applyQuaternion( quaternion );
			_this.object.up.applyQuaternion( quaternion );

			_rotateEnd.applyQuaternion( quaternion );

			if ( _this.staticMoving ) {

				_rotateStart.copy( _rotateEnd );

			} else {

				quaternion.setFromAxisAngle( axis, angle * ( _this.dynamicDampingFactor - 1.0 ) );
				_rotateStart.applyQuaternion( quaternion );

			}

		}

	};

	this.zoomCamera = function () {

		if ( _state === STATE.TOUCH_ZOOM ) {

			var factor = _touchZoomDistanceStart / _touchZoomDistanceEnd;
			_touchZoomDistanceStart = _touchZoomDistanceEnd;
			_eye.multiplyScalar( factor );

		} else {

			var factor = 1.0 + ( _zoomEnd.y - _zoomStart.y ) * _this.zoomSpeed;

			if ( factor !== 1.0 && factor > 0.0 ) {

				_eye.multiplyScalar( factor );

				if ( _this.staticMoving ) {

					_zoomStart.copy( _zoomEnd );

				} else {

					_zoomStart.y += ( _zoomEnd.y - _zoomStart.y ) * this.dynamicDampingFactor;

				}

			}

		}

	};

	this.panCamera = function () {

		var mouseChange = _panEnd.clone().sub( _panStart );

		if ( mouseChange.lengthSq() ) {

			mouseChange.multiplyScalar( _eye.length() * _this.panSpeed );

			var pan = _eye.clone().cross( _this.object.up ).setLength( mouseChange.x );
			pan.add( _this.object.up.clone().setLength( mouseChange.y ) );

			_this.object.position.add( pan );
			_this.target.add( pan );

			if ( _this.staticMoving ) {

				_panStart = _panEnd;

			} else {

				_panStart.add( mouseChange.subVectors( _panEnd, _panStart ).multiplyScalar( _this.dynamicDampingFactor ) );

			}

		}

	};

	this.checkDistances = function () {

		if ( !_this.noZoom || !_this.noPan ) {

			if ( _this.object.position.lengthSq() > _this.maxDistance * _this.maxDistance ) {

				_this.object.position.setLength( _this.maxDistance );

			}

			if ( _eye.lengthSq() < _this.minDistance * _this.minDistance ) {

				_this.object.position.addVectors( _this.target, _eye.setLength( _this.minDistance ) );

			}

		}

	};

	this.update = function () {

		_eye.subVectors( _this.object.position, _this.target );

		if ( !_this.noRotate ) {

			_this.rotateCamera();

		}

		if ( !_this.noZoom ) {

			_this.zoomCamera();

		}

		if ( !_this.noPan ) {

			_this.panCamera();

		}

		_this.object.position.addVectors( _this.target, _eye );

		_this.checkDistances();

		_this.object.lookAt( _this.target );

		if ( lastPosition.distanceToSquared( _this.object.position ) > 0 ) {

			_this.dispatchEvent( changeEvent );

			lastPosition.copy( _this.object.position );

		}

	};

	// listeners

	function keydown( event ) {

		if ( _this.enabled === false ) return;

		window.removeEventListener( 'keydown', keydown );

		_prevState = _state;

		if ( _state !== STATE.NONE ) {

			return;

		} else if ( event.keyCode === _this.keys[ STATE.ROTATE ] && !_this.noRotate ) {

			_state = STATE.ROTATE;

		} else if ( event.keyCode === _this.keys[ STATE.ZOOM ] && !_this.noZoom ) {

			_state = STATE.ZOOM;

		} else if ( event.keyCode === _this.keys[ STATE.PAN ] && !_this.noPan ) {

			_state = STATE.PAN;

		}

	}

	function keyup( event ) {

		if ( _this.enabled === false ) return;

		_state = _prevState;

		window.addEventListener( 'keydown', keydown, false );

	}

	function mousedown( event ) {

		if ( _this.enabled === false ) return;

		event.preventDefault();
		event.stopPropagation();

		if ( _state === STATE.NONE ) {

			_state = event.button;

		}

		if ( _state === STATE.ROTATE && !_this.noRotate ) {

			_rotateStart = _rotateEnd = _this.getMouseProjectionOnBall( event.clientX, event.clientY );

		} else if ( _state === STATE.ZOOM && !_this.noZoom ) {

			_zoomStart = _zoomEnd = _this.getMouseOnScreen( event.clientX, event.clientY );

		} else if ( _state === STATE.PAN && !_this.noPan ) {

			_panStart = _panEnd = _this.getMouseOnScreen( event.clientX, event.clientY );

		}

		document.addEventListener( 'mousemove', mousemove, false );
		document.addEventListener( 'mouseup', mouseup, false );

	}

	function mousemove( event ) {

		if ( _this.enabled === false ) return;

		event.preventDefault();
		event.stopPropagation();

		if ( _state === STATE.ROTATE && !_this.noRotate ) {

			_rotateEnd = _this.getMouseProjectionOnBall( event.clientX, event.clientY );

		} else if ( _state === STATE.ZOOM && !_this.noZoom ) {

			_zoomEnd = _this.getMouseOnScreen( event.clientX, event.clientY );

		} else if ( _state === STATE.PAN && !_this.noPan ) {

			_panEnd = _this.getMouseOnScreen( event.clientX, event.clientY );

		}

	}

	function mouseup( event ) {

		if ( _this.enabled === false ) return;

		event.preventDefault();
		event.stopPropagation();

		_state = STATE.NONE;

		document.removeEventListener( 'mousemove', mousemove );
		document.removeEventListener( 'mouseup', mouseup );

	}

	function mousewheel( event ) {

		if ( _this.enabled === false ) return;

		event.preventDefault();
		event.stopPropagation();

		var delta = 0;

		if ( event.wheelDelta ) { // WebKit / Opera / Explorer 9

			delta = event.wheelDelta / 40;

		} else if ( event.detail ) { // Firefox

			delta = - event.detail / 3;

		}

		_zoomStart.y += ( 1 / delta ) * 0.05;

	}

	function touchstart( event ) {

		if ( _this.enabled === false ) return;

		switch ( event.touches.length ) {

			case 1:
				_state = STATE.TOUCH_ROTATE;
				_rotateStart = _rotateEnd = _this.getMouseProjectionOnBall( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

			case 2:
				_state = STATE.TOUCH_ZOOM;
				var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
				var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
				_touchZoomDistanceEnd = _touchZoomDistanceStart = Math.sqrt( dx * dx + dy * dy );
				break;

			case 3:
				_state = STATE.TOUCH_PAN;
				_panStart = _panEnd = _this.getMouseOnScreen( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

			default:
				_state = STATE.NONE;

		}

	}

	function touchmove( event ) {

		if ( _this.enabled === false ) return;

		event.preventDefault();
		event.stopPropagation();

		switch ( event.touches.length ) {

			case 1:
				_rotateEnd = _this.getMouseProjectionOnBall( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

			case 2:
				var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
				var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
				_touchZoomDistanceEnd = Math.sqrt( dx * dx + dy * dy )
				break;

			case 3:
				_panEnd = _this.getMouseOnScreen( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

			default:
				_state = STATE.NONE;

		}

	}

	function touchend( event ) {

		if ( _this.enabled === false ) return;

		switch ( event.touches.length ) {

			case 1:
				_rotateStart = _rotateEnd = _this.getMouseProjectionOnBall( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

			case 2:
				_touchZoomDistanceStart = _touchZoomDistanceEnd = 0;
				break;

			case 3:
				_panStart = _panEnd = _this.getMouseOnScreen( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY );
				break;

		}

		_state = STATE.NONE;

	}

	this.domElement.addEventListener( 'contextmenu', function ( event ) { event.preventDefault(); }, false );

	this.domElement.addEventListener( 'mousedown', mousedown, false );

	this.domElement.addEventListener( 'mousewheel', mousewheel, false );
	this.domElement.addEventListener( 'DOMMouseScroll', mousewheel, false ); // firefox

	this.domElement.addEventListener( 'touchstart', touchstart, false );
	this.domElement.addEventListener( 'touchend', touchend, false );
	this.domElement.addEventListener( 'touchmove', touchmove, false );

	window.addEventListener( 'keydown', keydown, false );
	window.addEventListener( 'keyup', keyup, false );

	this.handleResize();

};

//#---------- end of file: "trackballcontrols.js"

//#---------- file: "three_view.js"

/*
  ECMAScript rendering in 3D of the solution to the ramp problem
  described in ./index.html
  

  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

/*global threeViewUpdate faceMesh document try_to_solve THREE threeViewUpdate*/

// Requires: ./three.js  and  ./solve_quartic.js

threeViewUpdate();  // at page load

function threeViewUpdate()
{
    var cont = document.getElementById( 'three-view' );

    var FLAT_VIEW = false;  // `true` used only to prepare the graph describing the problem in the horizontal plane.
    
    // --- fetch problem parameters

    var form = document.forms[ 'problem-input' ];

    var   depth = getFormFloat( form, 'red-depth' )
    , rampslope = getFormFloat( form, 'blue-slope' )
    , problem = {
        depth         : depth
        , horizlength : depth * 100 / rampslope  // in the horizontal plane
        , qq          : depth / Math.tan( Math.PI / 180 * getFormFloat( form, 'red-slope' ) )  // in the horizontal plane
        , width       : getFormFloat( form, 'blue-width' )  // just happens to have the same value as depth
    };

    var lp = threeViewUpdate.lastproblem;
    if (lp)
    {
        var different = false;
        for (var k in problem) { if (problem.hasOwnProperty( k )) {
            if (problem[ k ] !== lp[ k ])
            {
                different = true;
                break;
            }
        }}
        if (!different)
            return;
    }
    threeViewUpdate.lastproblem = problem;

    function getFormFloat( form, name )
    {
        return parseFloat( form.elements[ name ].value );
    }

    if (!threeViewUpdate.listening)
    {
        for (var k in form.elements)
        {
            var elt = form.elements[ k ];
            if (/^input$/i.test( elt.tagName ))
                elt.addEventListener( 'change', threeViewUpdate );  // on change
        }
        threeViewUpdate.listening = true;
    }

    // --- try to solve problem (== solve a quartic equation)

    var solution = try_to_solve( problem );

    // --- setup 3d view

    var rect = cont.getBoundingClientRect()
    ,   W    = rect.width | 0
    ,   H    = rect.height | 0
    , no_canvas_msg = '3D view not possible: Your browser does not support HTML5 Canvas. Think of trying browsers like Chrome, Safari or Firefox.'
    , renderer
    ;

    try {
        renderer = 'renderer' in threeViewUpdate  ?  threeViewUpdate.renderer  :  (threeViewUpdate.renderer = new THREE.CanvasRenderer);
    } 
    catch (e) 
    {
        alert( no_canvas_msg );
        renderer = threeViewUpdate.renderer = false;
    }
    
    if (!renderer)
    {
        cont.style.paddingTop = '30%';
        cont.style.boxSizing  = 'border-box';
        cont.innerHTML = no_canvas_msg;
    }
    else
    {
        renderer.setSize( W, H );
        
        if (renderer.domElement.parentNode !== cont)
            cont.appendChild( renderer.domElement );
        
        var scene = threeViewUpdate.scene  ||  (threeViewUpdate.scene = new THREE.Scene);
        
        var camera = threeViewUpdate.camera  ||  (threeViewUpdate.camera = new THREE.PerspectiveCamera(
            35,         // Field of view
            W / H,  // Aspect ratio
                .1,         // Near
            10000       // Far
        ));
    }

    // --- fetch solution (if any)
    
    var ar_deg = NaN
    , length   = NaN
    , gsslen   = NaN
    ;

    if (solution)
    {
        ar_deg = solution.angle_degrees;
        length = problem.horizlength * Math.sqrt( 1 + rampslope*rampslope / 1e4 );
        gsslen = solution.yy;
    }
    
    // --- populate result values in the article

    var out_elts = document.forms.result.elements;
    out_elts[ 'blue-length' ].value = '' + length;
    out_elts[ 'blue-angle'  ].value = '' + ar_deg;
    out_elts[ 'green-small-side-length' ].value = '' + gsslen;

    // xxx DOM message    if (!solution)
    

    // --- populate 3d view

    if (renderer)
    {
        if (FLAT_VIEW)
            depth = 0;

        while (scene.children.length)
            scene.remove( scene.children.slice( -1 )[ 0 ] );

        if (solution)
        {
            var   qq = problem.qq
            ,  width = problem.width
            , horizlength = problem.horizlength

            , xx = solution.xx
            , yy = solution.yy
            , zz = solution.zz
            , ar = solution.angle_radians

            , cut = Math.max( qq * 4, horizlength * Math.cos( ar ) * 1.5 )

            ;

            // Make sure we have all numbers
            [ cut, qq, width, horizlength, xx, yy, zz, ar ].forEach( function (x) { x.toPrecision.call.a; } );
            
            // Create a few faces

            scene.add( faceMesh( [ 0, { dx: +cut }, { dy: +qq, dz: -depth }, { x: qq } ]  ,  { color : 0xff0000 } ) );

            scene.add( faceMesh( [ 0, { x: qq, y: qq, z: -depth }, { y: cut }, { x: 0, z: 0 } ]  ,  { color : 0xff0000 } ) );

            var      cos_ar = Math.cos( ar )
            ,        sin_ar = Math.sin( ar )
            ,     yy_cos_ar = yy * cos_ar
            ,     yy_sin_ar = yy * sin_ar
            , horizlength_cos_ar = horizlength * cos_ar
            , horizlength_sin_ar = horizlength * sin_ar
            ,  width_cos_ar = width * cos_ar
            ,  width_sin_ar = width * sin_ar

            , corner_bottom = { x : yy_cos_ar  ,  y : yy_sin_ar }
            ,   left_point  = { x  : 0, y : zz, z : 0 }
            ;
            
            // Bottom

            scene.add( faceMesh(
                [ 
                    { x : qq, y : qq, z : -depth }
                    , { x : cut }
                    , { y : cut }
                    , { x : qq }
                ]
                , { color : 0x777777 }
            ) );

            // Ramp

            scene.add( faceMesh( 
                [ 
                    corner_bottom
                    , { dx : horizlength_cos_ar, dy : horizlength_sin_ar, dz : -depth }
                    , { dx : -width_sin_ar, dy : width_cos_ar }
                    , left_point
                ]
                ,  { color : 0x0000ff } 
            ) );

            // Fill the little triangular hole between borders and ramp

            scene.add( faceMesh(
                [
                    0
                    , corner_bottom
                    , left_point
                ]
                , { color : 0x22aa00 }
            ) )
        }
        
        // --- camera

        var to = new THREE.Vector3( cut/3, cut/3, 0 );

        if (!threeViewUpdate.cameraInitialized)
        {
            threeViewUpdate.cameraInitialized = true;
            
            camera.lookAt( to );
            
            camera.position.set( cut * 1.5, cut * 1.5, cut / 1.5 );
            camera.up.set( 0, 0, 1 );
        }
        
        scene.add( camera );
        
        // --- light
        
        var light = new THREE.PointLight( 0xFFFFFF );
        light.position.set( cut/2, cut/5, 100 );
        scene.add( light );   

        // --- render 3d view
        
        render();

        // --- animate 3d view: rotation through mouse "drag and drop"

        if (!threeViewUpdate.controls)
        {
            var controls = threeViewUpdate.controls = new THREE.TrackballControls( camera, renderer.domElement );
            controls.target.set( to.x, to.y, 0 )
            
            controls.rotateSpeed = 1.0;
            controls.zoomSpeed = 1.2;
            controls.panSpeed = 0.8;
            
            controls.noZoom = false;
            controls.noPan = false;
            
            controls.staticMoving = true;
            controls.dynamicDampingFactor = 0.3;
            
            controls.keys = [ 65, 83, 68 ];
            
            controls.addEventListener( 'change', render );

            animate();

        }
    }

    // --- Details (function declarations)

    function render()
    {
        renderer.render( scene, camera );
    }


    function animate() 
    {
	requestAnimationFrame( animate );
	controls.update();
    }
}

function faceMesh( arr, opt )
// Helper to create a simple mesh.
// 
// Default rotation: -90 degrees on both x & z axes, for better mouse
// rotation usability in the above use case (TrackballControls).
{
    var cfg = opt  ||  {};

    var color = cfg.color  ||  0xFF0000
    ,   rot_x = cfg.rot_x != null  ?  cfg.rot_x  :  0
    ,   rot_y = cfg.rot_y != null  ?  cfg.rot_y  :  0
    ,   rot_z = cfg.rot_z != null  ?  cfg.rot_z  :  0
    ;
    
    var geometry = new THREE.Geometry();
    
    for (var x = 0, y = 0, z = 0, n = arr.length, 
         i = 0;
         i < n; 
         i++ )
    {
        var spec = arr[ i ];
        if (spec == 0)
        {
            x = y = z = 0;
        }
        else
        {
            if (spec.dx)   
                x += spec.dx;
            else if (spec.x != null)
                x = spec.x;
            
            if (spec.dy) 
                y += spec.dy;
            else if (spec.y != null)
                y = spec.y;
            
            if (spec.dz) 
                z += spec.dz;
            else if (spec.z != null)
                z = spec.z;
        }
        
        geometry.vertices.push( new THREE.Vector3( x, y, z ) );
    }

    if (n == 3)
    {
        geometry.faces.push( new THREE.Face3( 0, 1, 2 ) );
        geometry.faces.push( new THREE.Face3( 2, 1, 0 ) );
    }
    else if (n == 4)
    {
        geometry.faces.push( new THREE.Face4( 0, 1, 2, 3 ) );
        geometry.faces.push( new THREE.Face4( 3, 2, 1, 0 ) );
    }
    else
    {
        throw new Error( 'n ' + n + ' not supported.' );
    }

    geometry.computeFaceNormals();
    
    var ret = new THREE.Mesh(
        geometry,
        new THREE.MeshLambertMaterial( { color: color } )
    );

    ret.rotation.x = rot_x;
    ret.rotation.y = rot_y;
    ret.rotation.z = rot_z;

    return ret;
}

//#---------- end of file: "three_view.js"

