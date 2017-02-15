var express = require('express')

var bodyParser = require('body-parser');

var app = express();


app.use(bodyParser.json());
app.use(bodyParser.urlencoded({
	extended: true
}));


app.get('/', function(req,res){
	res.sendFile(__dirname+ "/" + "form.html")
})

app.post('/', function(req,res){
	var first = parseInt(req.body.a,10);
	var second = parseInt(req.body.b,10);
	var answer = first +second;
	console.log("first is = " ,parseInt(first));
	res.send("answer is " + answer);


})


var server = app.listen(8081, function(){
	var host = server.address().address
	var port = server.address().port

	console.log("running at http://%s:%s", host,port)
})

