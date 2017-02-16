$(document).foundation()

// global variables
var testVar = "global";
var outBoxContent = "This is output box text.";


// functions
function uploadSourceGrid(){
   return("Source grid is uploaded.");
}


// document ready event
$(document).ready(function(){

  // sample hello world alert
  $("#hello").click(function() {
    alert("Hello, world!");
  });

  // populate the outbox
  $('#gridStats').click(function(){
    $('#outbox').val(outBoxContent);
  });

  // transfer button click event
  $('#transfer').bind('click', function( event ){
     alert('Performing transfer action.');
  });

  // srcUploader button change event
  $('#srcGridFileUpload').bind('change', function( event ){
   $('#outbox').val( uploadSourceGrid() );
  });

  // dummy tests
  //$("a").css("background-color", "blue");

});
