
// the first row of the 2d array needs to contain the values of the object keys

// into an object of named vectors
function convert_2d_array_into_df_object(arr_) {

  // create a dataframe-like object df={"column_header1":[[],[]]}
  // each key/column contains a vector (still a 2d-array)
  let df = {};
  for (let x = 0; x < arr_[0].length; x++) {
    
    let column_name = arr_[0][x];
    if (!df[column_name]) {
      df[column_name] = [];
    }
    
    for (let y = 1; y < arr_.length; y++) {
      
      let value = parseFloat(arr_[y][x]);
      df[column_name].push([ value ]);
    }

  }
  return df;
}; // closing fn


function convert_2d_array_to_arr_of_objs(arr_) {
  
  let output_arr_ = [];

  for (let y = 1; y < arr_.length; y++) {
    let obj = {};
    for (let x = 0; x < arr_[0].length; x++) {
      // let key = string_arr[arr_[0][x]];
      let key = arr_[0][x];
      obj[key] = parseFloat(arr_[y][x]);
    } // closing x-loop
    output_arr_.push(obj);
  } // closing y-loop
  return output_arr_;
}; // closing fn