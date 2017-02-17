// globals
var express =   require("express");
var app     =   express();
var path = require('path');
var formidable = require('formidable');
var fs = require('fs');

// other globals
var verb;
var serverStartupDate = new Date();


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


app.use(express.static(path.join(__dirname, '../public')));

app.get('/', function(req, res){
  res.sendFile(path.join(__dirname, '../public/index.html'));
});

app.post('/upload', function(req, res){
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
    fs.rename(file.path, path.join(form.uploadDir, "source.stl"));
    console.log("File received: %s", path.join(form.uploadDir, file.name));
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

app.listen(3000,function(){
    console.log("Working on port 3000");
});
