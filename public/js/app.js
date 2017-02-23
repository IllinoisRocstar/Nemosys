$(document).foundation()

// global variables
// nothing so far


// document ready event
$(document).ready(function(){


  // sample hello world alert
  $("#hello").click(function() {
    alert("Hello, world!");
  });


  // source grid statistics button
  $('#gridStats').click(function(){
    $('#outbox').val("Submitting source grid statistics request");
    // create CMLHttpRequest object, ActiveX object if old IE5, IE6
    var xhttp;
    if (window.XMLHttpRequest) {
      xhttp = new XMLHttpRequest();
    } else {
      // code for IE6, IE5
      xhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    // preparing and sending a get request
    $('#outbox').val("Sending request to server.\n");
    xhttp.onreadystatechange = function() {
      if (this.readyState == 4 && this.status == 200) {
        $('#outbox').val('Source Grid Statistics \n'+
                         '------------------------------ \n' + 
                         this.responseText);
      }
    };
    xhttp.open("GET", "srcGrdStats", true);
    xhttp.send();
  });

  // source grid statistics button
  $('#gridSlnNames').click(function(){
    $('#outbox').val("Submitting source grid solution names request");
    // create CMLHttpRequest object, ActiveX object if old IE5, IE6
    var xhttp;
    if (window.XMLHttpRequest) {
      xhttp = new XMLHttpRequest();
    } else {
      // code for IE6, IE5
      xhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    // preparing and sending a get request
    $('#outbox').val("Sending request to server.\n");
    xhttp.onreadystatechange = function() {
      if (this.readyState == 4 && this.status == 200) {
        $('#outbox').val('Source Grid Solution Names \n'+
                         '------------------------------------------ \n' + 
                         this.responseText);
      }
    };
    xhttp.open("GET", "srcGrdSlnNames", true);
    xhttp.send();
  });


  // transfer solution button click event
  $('#transferSln').bind('click', function( event ){
    $('#outbox').val("Submitting solution transfer request");
    // create CMLHttpRequest object, ActiveX object if old IE5, IE6
    var xhttp;
    if (window.XMLHttpRequest) {
      xhttp = new XMLHttpRequest();
    } else {
      // code for IE6, IE5
      xhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    // preparing and sending a get request
    $('#outbox').val("Sending request to server.\n");
    xhttp.onreadystatechange = function() {
      if (this.readyState == 4 && this.status == 200) {
        $('#outbox').val('Solution Transfer Log \n'+
                         '----------------------------- \n' + 
                         this.responseText);
      }
    };
    xhttp.open("GET", "slnTransfer", true);
    xhttp.send();
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
       url: '/uploadSrc',
       type: 'POST',
       data: formData,
       processData: false,
       contentType: false,
       success: function(data){
	   console.log('upload successful!\n' + data);
           initialize();
           $('#outbox').val('Source uploaded successfully!\n');
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


  // trgUploader button change event
  $('#trgGridFileUpload').bind('change', function(){    

   // zeroing progress bar
   $('#trgUploadProgBar').text('0%');
   $('#trgUploadProgBar').width('0%');
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
       url: '/uploadTrg',
       type: 'POST',
       data: formData,
       processData: false,
       contentType: false,
       success: function(data){
	   console.log('Target uploaded successful!\n' + data);
           initialize();
           $('#outbox').val('Target uploaded successfully!\n');
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
	     $('#trgUploadProgBar').text(percentComplete + '%');
	     $('#trgUploadProgBar').width(percentComplete + '%');
	     // once the upload reaches 100%, set the progress bar text to done
	     if (percentComplete === 100) {
	       $('#trgUploadProgBar').html('Done');
	     }
	   }
	 }, false);
         return xhr;
       }
     });
   }    

  });

  // grid download button
  $("#trgDownload").click(function() {
      window.location = './uploads/fluid_06.100000_0000.cgns';
  });

  // dummy tests
  //$("a").css("background-color", "blue");

});
