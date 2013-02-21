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

