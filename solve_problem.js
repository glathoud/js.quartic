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
