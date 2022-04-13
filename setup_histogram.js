(function() {
  
  window.addEventListener('load', function() {

    let arr_x = [];
    for (let i = 0; i < 10000; i++) {
      arr_x.push(chisq(5));  
    }
    
    let hist_x_minimum = new Histogram({'arr':arr_x});
    container_histogram_minimum.appendChild(hist_x_minimum.RETURN_CANVAS());
    
    let hist_x_options = new Histogram({'arr':arr_x,'number_of_bins':23});
    container_histogram_options.appendChild(hist_x_options.RETURN_CANVAS());
    

  }); // CLOSING window.onload

})(); // CLOSING anon