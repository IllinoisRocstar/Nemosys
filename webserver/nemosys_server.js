// globals

var express =   require("express");
var app     =   express();
var path = require('path');
var formidable = require('formidable');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));
var ip;



const spawnMe = require('child_process').spawn;   // using a new shell
const execMe = require('child_process').execFile; // using current shell


// other globals
var verb;
var serverStartupDate = new Date();
var g2gCmd = '/home/msafdari/scratch/scratch/Nemosys/build/bin/grid2gridTransfer';
var scratchFldr = '/home/msafdari/scratch/scratch/WebNemosys/webserver/scratch';


// passed command line arguments
var args = process.argv.slice(2);
console.log("Welcome to NEMoSys Webserver Backend(v1.0)");
console.log("Starting server in : %s", serverStartupDate);


// process input switches
function processArgs(){
   //console.log("The server called with these switches : ", args);
   verb = argv['v'];
   // grid to grid conversion commands
   if (argv['g2g']) {
     g2gCmd =  argv['g2g'];
     console.log('grid2grid application : ' + g2gCmd);
   }
   // scratch folder location
   if (argv['scratch']) {
     scratchFldr =  argv['scratch'];
     console.log('scratch folder : ' + scratchFldr);
   }
}
processArgs();

// setting up the path
app.use(express.static(path.join(__dirname, '../public')));


// webpage requests
app.get('/', function(req, res){
  res.sendFile(path.join(__dirname, '../public/index.html'));
});





// upload source grid
app.post('/uploadSrc', function(req, res){
  // create an incoming form object
  var form = new formidable.IncomingForm();
  // specify that we want to allow the user to upload multiple files in a single request
  form.multiples = true;
  // store all uploads in the /uploads directory
  form.uploadDir = path.join(__dirname, '../public/uploads');
  // every time a file has been uploaded successfully,
  // rename it to it's orignal name


  form.on('file', function(field, file) {
    fs.rename(file.path, path.join(form.uploadDir, file.name));
    //fs.rename(file.path, path.join(form.uploadDir,  "fluid_04.100000_0000.cgns"));
    console.log("File received: %s", path.join(form.uploadDir, file.name));
    // skin the CGNS mesh
    execMe(g2gCmd, ['--cgns2stl','../../public/uploads/fluid_04.100000_0000.cgns', 
                    '../../public/uploads/', 'src'], {cwd: scratchFldr}, 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`${stderr}`);
	  console.log(`Error: ${error}`);
	} 
	else
	{
	  console.log(`${stdout}`);
          res.end('success');
	}
     });

  });
  // log any errors that occur
  form.on('error', function(err) {
    console.log('An error has occured: \n' + err);
  });
  // once all the files have been uploaded, send a response to the client
  //form.on('end', function() {
    //res.end('success');
  //});
  // parse the incoming request containing the form data
  form.parse(req);
});

// upload target grid
app.post('/uploadTrg', function(req, res){
  // create an incoming form object
  var form = new formidable.IncomingForm();
  // specify that we want to allow the user to upload multiple files in a single request
  form.multiples = true;
  // store all uploads in the /uploads directory
  form.uploadDir = path.join(__dirname, '../public/uploads');
  // every time a file has been uploaded successfully,
  // rename it to it's orignal name
  form.on('file', function(field, file) {
    //fs.rename(file.path, path.join(form.uploadDir, file.name));
    fs.rename(file.path, path.join(form.uploadDir, "fluid_06.100000_0001.cgns"));
    console.log("File received: %s", path.join(form.uploadDir, file.name));
    // skin the CGNS mesh
    execMe(g2gCmd,  ['--cgns2stl','../../public/uploads/fluid_06.100000_0001.cgns', 
                    '../../public/uploads/', 'trg'],  {cwd: scratchFldr}, 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`${stderr}`);
	  console.log(`Error: ${error}`);
	} 
	else
	{
	  console.log(`${stdout}`);
          res.end('success');
	}
     });

  });
  // log any errors that occur
  form.on('error', function(err) {
    console.log('An error has occured: \n' + err);
  });
  // once all the files have been uploaded, send a response to the client
  //form.on('end', function() {
    //res.end('success');
  //});
  // parse the incoming request containing the form data
  form.parse(req);
});


// statistics of source grid
app.get('/srcGrdStats', function(req, res){
  console.log("Source grid statistics request received.\n");

     execMe(g2gCmd, ['--statCGNS','../../public/uploads/fluid_04.100000_0000.cgns'],  {cwd: scratchFldr}, 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`${stderr}`);
	  console.log(`Error:\n ${error}`);
	  res.send('Error :\n ${error}')
	} 
	else
	{
	  console.log(`${stdout}`);
	  res.send(`${stdout}`);
	}
     });
});


// solution names of source grid
app.get('/srcGrdSlnNames', function(req, res){
  console.log("Source grid solution name request received.\n");

     execMe(g2gCmd, ['--listSln','../../public/uploads/fluid_04.100000_0000.cgns'],  {cwd: scratchFldr}, 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`${stderr}`);
	  console.log(`Error:\n ${error}`);
	  res.send('Error :\n ${error}')
	} 
	else
	{
	  console.log(`${stdout}`);
	  res.send(`${stdout}`);
	}
     });
});


// solution transfer from source grid to target
app.get('/slnTransfer', function(req, res){
  console.log("Solution transfer request received.\n");

     execMe(g2gCmd, ['--transCGNS', '../../public/uploads/fluid_04.100000_0000.cgns',
                     '../../public/uploads/fluid_06.100000_0000.cgns', 
                     '../../public/uploads/target_with_sln.cgns'],  {cwd: scratchFldr}, 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`${stderr}`);
	  console.log(`Error:\n ${error}`);
	  res.send('Error :\n ${error}')
	} 
	else
	{
	  console.log(`${stdout}`);
	  res.send(`${stdout}`);
	}
     });
});


// running the webserver
app.listen(3000,function(){
    console.log("Working on port 3000");
});
