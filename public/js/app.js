$(document).foundation()

// global variables
var testVar = "global";
var outBoxContent = "This is output box text.";


// functions
function uploadSourceGrid(){
   return("Source grid is uploaded successfully.");
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
  $('#srcGridFileUpload').bind('change', function(){    

   // zeroing progress bar
   $('#srcUploadProgBar').text('0%');
   $('#srcUploadProgBar').width('0%');
   var files = $(this).get(0).files;
   if (files.length > 0){
     // One or more files selected, process the file upload
     // create a FormData object which will be sent as the data payload in the
     // AJAX request
     var formData = new FormData();
     // loop through all the selected files
     for (var i = 0; i < files.length; i++) {
       var file = files[i];
       // add the files to formData object for the data payload
       formData.append('uploads[]', file, file.name);
     }
     $.ajax({
       url: '/upload',
       type: 'POST',
       data: formData,
       processData: false,
       contentType: false,
       success: function(data){
	   console.log('upload successful!\n' + data);
           $('#outbox').val('Source uploaded successfully!\n' + data);
       },
       xhr: function() {
	 // create an XMLHttpRequest
	 var xhr = new XMLHttpRequest();
	 // listen to the 'progress' event
	 xhr.upload.addEventListener('progress', function(evt) {
	   if (evt.lengthComputable) {
	     // calculate the percentage of upload completed
	     var percentComplete = evt.loaded / evt.total;
	     percentComplete = parseInt(percentComplete * 100);
	     // update the Bootstrap progress bar with the new percentage
	     $('#srcUploadProgBar').text(percentComplete + '%');
	     $('#srcUploadProgBar').width(percentComplete + '%');
	     // once the upload reaches 100%, set the progress bar text to done
	     if (percentComplete === 100) {
	       $('#srcUploadProgBar').html('Done');
	     }
	   }
	 }, false);
         return xhr;
       }
     });
   }    

  });



  // dummy tests
  //$("a").css("background-color", "blue");

});
