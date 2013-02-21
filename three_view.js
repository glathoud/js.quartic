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
