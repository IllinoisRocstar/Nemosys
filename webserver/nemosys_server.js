// globals
var express =   require("express");
var app     =   express();
var path = require('path');
var formidable = require('formidable');
var fs = require('fs');


const spawnMe = require('child_process').spawn;   // using a new shell
const execMe = require('child_process').execFile; // using current shell

//const ls = spawnMe('ls', ['-lh', '/usr']);


// other globals
var verb;
var serverStartupDate = new Date();
var g2gCmd = '/home/msafdari/scratch/scratch/Nemosys/build/bin/grid2gridTransfer';


// passed command line arguments
var args = process.argv.slice(2);
console.log("Welcome to NEMoSys Webserver Backend(v1.0)");
console.log("Starting server in : %s", serverStartupDate);
console.log("The server called with these switches : ", args);


// process input switches
process.argv.forEach(function (val, index, array) {
  if (val=='-v') {    
    console.log("Working in the verbose mode.");
    verb = true;
  }
});


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
    //fs.rename(file.path, path.join(form.uploadDir, file.name));
    fs.rename(file.path, path.join(form.uploadDir, "fluid_04.100000_0000.cgns"));
    console.log("File received: %s", path.join(form.uploadDir, file.name));
    // skin the CGNS mesh
    execMe(g2gCmd, ['--cgns2stl','../public/uploads/fluid_04.100000_0000.cgns', '../public/uploads/'], 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`stdout: ${stderr}`);
	  console.log(`error: ${error}`);
	} 
	else
	{
	  console.log(`stdout: ${stdout}`);
	}
     });

  });
  // log any errors that occur
  form.on('error', function(err) {
    console.log('An error has occured: \n' + err);
  });
  // once all the files have been uploaded, send a response to the client
  form.on('end', function() {
    res.end('success');
  });
  // parse the incoming request containing the form data
  form.parse(req);
});


// statistics of source grid
app.get('/srcGrdStats', function(req, res){
  console.log("Source grid statistics request received.\n");

  var g2gtNemo = execMe(g2gCmd, ['--statCGNS','../public/uploads/fluid_04.100000_0000.cgns'], 
     function callback(error, stdout, stderr) {
	if (error){
	  console.log(`stdout: ${stderr}`);
	  console.log(`error: ${error}`);
	  res.send('error : ${error}')
	} 
	else
	{
	  console.log(`stdout: ${stdout}`);
	  res.send(`stdout: ${stdout}`);
	}
     });
});


// running the webserver
app.listen(3000,function(){
    console.log("Working on port 3000");
});
