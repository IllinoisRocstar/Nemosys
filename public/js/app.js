$(document).foundation()

// global variables
// nothing so far
var sourceFormatted;
var targetFormatted;
var targetSln;
// document ready event
$(document).ready(function(){


  // sample hello world alert
  $("#hello").click(function() {
    alert("Hello, world!");
  });

  //WK adding plot button
  // Generating Histogram via running nemosys backend synchronyously
  $('#plot').click(function(){
    if(typeof sourceFormatted == "undefined"){
    //alert("souce format = " + sourceFormatted)
    $('#outbox').val("No Source Grid selected");
    //$('#collapse').foundation('toggle')
    } else{
    $('#outbox').val("Submitting histogram.json generation request");
 
    // create CMLHttpRequest object, ActiveX object if old IE5, IE6
    var xhttp;
    if (window.XMLHttpRequest) {
      xhttp = new XMLHttpRequest();
    } else {
      // code for IE6, IE5
      xhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    
    // callback option = false
    var params = 'formatName='+sourceFormatted
    //alert(params)
    xhttp.open("GET", "plot"+"?"+params, false);
    xhttp.send();

    // Display Histogram
    var xhttpDisplayHist;
    var jsonfile;
    if (window.XMLHttpRequest) {
      xhttpDisplayHist = new XMLHttpRequest();
    } else {
      // code for IE6, IE5
      xhttpDisplayHist = new ActiveXObject("Microsoft.XMLHTTP");
    }
    // preparing and sending a get request
    $('#outbox').val("Sending request to server.\n");
    xhttpDisplayHist.onreadystatechange = function() {
      if (this.readyState == 4 && this.status == 200) {
        jsonfile = JSON.parse(this.response);
        $('#outbox').val('Histogram Statistics \n'+
                        '------------------------------ \n' + 
						"0.0  < Q <  0.1 :   " + jsonfile[0] + " elements\n" +
						"0.1  < Q <  0.2 :   " + jsonfile[1] + " elements\n" +
						"0.2  < Q <  0.3 :   " + jsonfile[2] + " elements\n" +
						"0.3  < Q <  0.4 :   " + jsonfile[3] + " elements\n" +
						"0.4  < Q <  0.5 :   " + jsonfile[4] + " elements\n" +
						"0.5  < Q <  0.6 :   " + jsonfile[5] + " elements\n" +
						"0.6  < Q <  0.7 :   " + jsonfile[6] + " elements\n" +
						"0.7  < Q <  0.8 :   " + jsonfile[7] + " elements\n" +
						"0.8  < Q <  0.9 :   " + jsonfile[8] + " elements\n" +
						"0.9  < Q <  1.0 :   " + jsonfile[9] + " elements\n"); 
      }
    };
    
    // callback option is true this time as we want to display the graph w/o refreshing
    xhttpDisplayHist.open("GET", "getHist", false);
    xhttpDisplayHist.send();


    // Google Chart code
    google.charts.load("current", {packages:['corechart']});
    google.charts.setOnLoadCallback(drawChart);

    var loc = window.location.pathname;

    function drawChart() {
      var data = google.visualization.arrayToDataTable([
        ["Element", "Number of Elements"],
        ["0.0  < Q <  0.1", parseInt(jsonfile[0])],
        ["0.1  < Q <  0.2", parseInt(jsonfile[1])],
        ["0.2  < Q <  0.3", parseInt(jsonfile[2])],
        ["0.3  < Q <  0.4", parseInt(jsonfile[3])],
        ["0.4  < Q <  0.5", parseInt(jsonfile[4])],
        ["0.5  < Q <  0.6", parseInt(jsonfile[5])],
        ["0.6  < Q <  0.7", parseInt(jsonfile[6])],
        ["0.7  < Q <  0.8", parseInt(jsonfile[7])],
        ["0.8  < Q <  0.9", parseInt(jsonfile[8])],
        ["0.9  < Q <  1.0", parseInt(jsonfile[9])],
      ]);
       
      // A lot of visualization option can be set here under option object
      var view = new google.visualization.DataView(data);
      var options = {
        title: "Histogram of element distribution",
        width: 1100,
        height: 400,
        hAxis : {
           textStyle : {
              fontSize:9
           }
        }
      };

      var chart = new google.visualization.ColumnChart(document.getElementById("chart_diz"));
      chart.draw(view, options);
  }
}
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
    var params = 'formatName='+sourceFormatted
    //alert(params)
    xhttp.open("GET", "srcGrdStats"+"?"+params, true);
    xhttp.send();
  });

  // source solution name button
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
    var params = 'formatName='+sourceFormatted
    xhttp.open("GET", "srcGrdSlnNames"+"?"+params, true);
    xhttp.send();
  });


  // transfer solution button click event
  $('#transferSln').bind('click', function( event ){
    if(typeof sourceFormatted == "undefined" || typeof targetFormatted == "undefined" ){
      $('#outbox').val("A pair of source and target grid first need to be selected");
      alert("A pair of source and target grid first need to be selected");
      return;
    }
    targetSln = './uploads/target_with_sln_'+targetFormatted;
    ip="_" +document.getElementById("ip").innerHTML + "_" ;
    targetSln = targetSln.replace(ip,'_');
    $('#outbox').val("Submitting solution transfer request");
    // checkbox value
    var params;
    if ($('#transErrChk').prop('checked')==true) {
       params = "checkSlnTransErr=true";
    } else {
       params = "checkSlnTransErr=false";
    }
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
    var src = 'formatSrcName='+sourceFormatted
    var trg = 'formatTrgName='+targetFormatted
    var ip = 'formatIPName='+ip
    xhttp.open("GET", "slnTransfer"+"?"+params+"&"+src+"&"+trg+"&"+ip, false);
    xhttp.send();
  });


  // srcUploader button change event
  $('#srcGridFileUpload').bind('change', function(){    
    var ip;
   // zeroing progress bar
   $('#srcUploadProgBar').text('0%');
   $('#srcUploadProgBar').width('0%');
   
   var files = $(this).get(0).files;
   if (files.length > 0){
     // One or more files selected, process the file upload
     // create a FormData object which will be sent as the data payload in the
     // AJAX request
     var formData = new FormData();
     
     ip=document.getElementById("ip").innerHTML;
     //window.alert("IP: " + ip);
     // loop through all the selected files
     for (var i = 0; i < files.length; i++) {
       var file = files[i];
       // add the files to formData object for the data payload
       var d = new Date();
       var newName = d.getDate() + "-" + (d.getMonth()+1) + "-" + d.getFullYear() + "_" + d.getTime() + "_" + ip + "_" + file.name;
       //window.alert(ip_name);
       formData.append('uploads[]', file,  newName);
     }
     sourceFormatted = newName;
     $.ajax({
       url: '/uploadSrc',
       type: 'POST',
       data: formData,
       processData: false,
       contentType: false,
       success: function(data){
	   console.log('upload successful!\n' + data);
           initialize(sourceFormatted,targetFormatted);
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

  // To close the panel when upload source grid is clicked
  $('.button').first().on('click',function(){
    $('#collapse').fadeOut(40);
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



     ip=document.getElementById("ip").innerHTML;
     // loop through all the selected files
     for (var i = 0; i < files.length; i++) {
       var file = files[i];
       // add the files to formData object for the data payload
       var d = new Date();
       var newName = d.getDate() + "-" + (d.getMonth()+1) + "-" + d.getFullYear() + "_" + d.getTime() + "_" + ip + "_" +file.name;
       formData.append('uploads[]', file, newName);
     }
     targetFormatted = newName;
     $.ajax({
       url: '/uploadTrg',
       type: 'POST',
       data: formData,
       processData: false,
       contentType: false,
       success: function(data){
	   console.log('Target uploaded successful!\n' + data);
           initialize(sourceFormatted,targetFormatted);
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
    if(typeof targetSln !== 'undefined'){
      $.ajax({
        url:targetSln,
        type:'HEAD',
        error: function(){
          alert("No solution file available.")
        },
        success : function(){
        window.location = targetSln;
        }

      });
    }else alert("Click 'Transfer Solution'");
  });

 // dummy tests
  //$("a").css("background-color", "blue");


});
