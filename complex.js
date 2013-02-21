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
