var stats;
var containerSrc, loaderSrc, cameraSrc, controlsSrc, sceneSrc, rendererSrc;
var containerTrg, loaderTrg, cameraTrg, controlsTrg, sceneTrg, rendererTrg;

initialize();
animate();


function initialize() {
  console.log("Starting initialize()");

  containerSrc = document.getElementById( 'webglContainerSrc' ) 
  containerTrg = document.getElementById( 'webglContainerTrg' ) 

  // create a scenes 
  sceneSrc = new THREE.Scene();
  sceneSrc.background = new THREE.Color( 0 );
  sceneTrg = new THREE.Scene();
  sceneTrg.background = new THREE.Color( 0 );

  // Create cameras which will be appened to the scene objects
  cameraSrc = new THREE.PerspectiveCamera( 60, 1.0, 0.01, 1e10 );
  cameraSrc.position.x = 0.01;
  cameraSrc.position.y = 0.01;
  cameraSrc.position.z = 0.01;
  cameraTrg = new THREE.PerspectiveCamera( 60, 1.0, 0.01, 1e10 );
  cameraTrg.position.x = 0.01;
  cameraTrg.position.y = 0.01;
  cameraTrg.position.z = 0.01;

  // Create controls object which allows for user control of the camera positions
  controlsSrc = new THREE.TrackballControls( cameraSrc , containerSrc );
  controlsSrc.rotateSpeed = 5.0;
  controlsSrc.zoomSpeed = 5;
  controlsSrc.panSpeed = 2;
  controlsSrc.noZoom = false;
  controlsSrc.noPan = false;
  controlsSrc.staticMoving = true;
  controlsSrc.dynamicDampingFactor = 0.3;
  controlsTrg = new THREE.TrackballControls( cameraTrg , containerTrg );
  controlsTrg.rotateSpeed = 5.0;
  controlsTrg.zoomSpeed = 5;
  controlsTrg.panSpeed = 2;
  controlsTrg.noZoom = false;
  controlsTrg.noPan = false;
  controlsTrg.staticMoving = true;
  controlsTrg.dynamicDampingFactor = 0.3;

  // append the cameras to the scene objects
  sceneSrc.add( cameraSrc );
  sceneTrg.add( cameraTrg );

  // adding lights
  var dirLightSrc = new THREE.DirectionalLight( 0xffffff ); // directional light
  dirLightSrc.position.set( 200, 200, 1000 ).normalize(); 
  cameraSrc.add( dirLightSrc );
  cameraSrc.add( dirLightSrc.target );
  var ambientLightSrc = new THREE.AmbientLight(0x0c0c0c); 
  cameraSrc.add(ambientLightSrc);
  var dirLightTrg = new THREE.DirectionalLight( 0xffffff ); // directional light
  dirLightTrg.position.set( 200, 200, 1000 ).normalize(); 
  cameraTrg.add( dirLightTrg );
  cameraTrg.add( dirLightTrg.target );
  var ambientLightTrg = new THREE.AmbientLight(0x0c0c0c); 
  cameraTrg.add(ambientLightTrg);

  // Create xyz axes at the origins
  var axesSrc = new THREE.AxisHelper(1); // length of axes displayed is argument 
  sceneSrc.add(axesSrc);
  var axesTrg = new THREE.AxisHelper(1); // length of axes displayed is argument 
  sceneTrg.add(axesTrg);

  // note: In threejs, in order to display an object you need to create a "mesh". A "mesh" 
  // consists of a "geometry" and a "material". 
  var materialSrc = new THREE.MeshLambertMaterial( 
      { 
	  color: 0xffffff, 
	  side: THREE.DoubleSide, 
	  wireframe: false 
      });
  var materialTrg = new THREE.MeshLambertMaterial( 
      { 
	  color: 0xfff111, 
	  side: THREE.DoubleSide, 
	  wireframe: false 
      });

  // note: the VTKLoader only works for a subset of VTK!!! 
  // Alternatively, the "STLLoader" can be used with the same syntax
  // var loaderSrc = new THREE.VTKLoader();
  loaderSrc = new THREE.STLLoader();
  loaderSrc.load( "uploads/source.stl", function ( geometry ) {
      geometry.center();
      geometry.computeVertexNormals();
      geometry.dynamic = true;
      geometry.verticesNeedUpdate = true;
      var meshSrc = new THREE.Mesh( geometry, materialSrc ); 
      meshSrc.position.set( 0, 0, 0 );
      meshSrc.scale.multiplyScalar( 1 );
      sceneSrc.add( meshSrc ); // append the mesh to the sceneSrc
  } );
  loaderTrg = new THREE.STLLoader();
  loaderTrg.load( "uploads/target.stl", function ( geometry ) {
      geometry.center();
      geometry.computeVertexNormals();
      geometry.dynamic = true;
      geometry.verticesNeedUpdate = true;
      var meshTrg = new THREE.Mesh( geometry, materialTrg ); 
      meshTrg.position.set( 0, 0, 0 );
      meshTrg.scale.multiplyScalar( 1 );
      sceneTrg.add( meshTrg ); // append the mesh to the sceneTrg
  } );

  // Create renderers
  rendererSrc = new THREE.WebGLRenderer( { antialias: true } );
  rendererSrc.setPixelRatio( window.devicePixelRatio );
  rendererSrc.setSize( 380, 380 );
  rendererTrg = new THREE.WebGLRenderer( { antialias: true } );
  rendererTrg.setPixelRatio( window.devicePixelRatio );
  rendererTrg.setSize( 380, 380 );

  // <----- in template, create a div with this ... repeat for additional containerSrc, ie 
  // one for the "input mesh" and one for the "target mesh"
  // rendererSrc.setSize( containerSrc.innerWidth, containerSrc.innerHeight ); 
  document.getElementById('containerSrc').appendChild( containerSrc );
  containerSrc.appendChild( rendererSrc.domElement );
  document.getElementById('containerTrg').appendChild( containerTrg );
  containerTrg.appendChild( rendererTrg.domElement );

  // display stats (not needed, debug only)
  stats = new Stats();
  containerSrc.appendChild( stats.dom );

  // watch for window resize - see function def below
  window.addEventListener( 'resize', onWindowResize, false );
}



function onWindowResize() {
  cameraSrc.aspect = window.innerWidth / window.innerHeight;
  cameraSrc.updateProjectionMatrix();
  rendererSrc.setSize( 379, 380 );
  controlsSrc.handleResize();
  cameraTrg.aspect = window.innerWidth / window.innerHeight;
  cameraTrg.updateProjectionMatrix();
  rendererTrg.setSize( 379, 380 );
  controlsTrg.handleResize();
}


// --- this actually renders the screen and attempts to maintain 
// a constant 60 fps refresh rate (ie, this is "called" 60 times per sec) ---------------
function animate() {
  rendererSrc.render( sceneSrc, cameraSrc );
  controlsSrc.update();
  rendererTrg.render( sceneTrg, cameraTrg );
  controlsTrg.update();
  stats.update();
  requestAnimationFrame( animate );
}

