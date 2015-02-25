//aligns
function align(seq1, seq2, seq3) {
  //initialize arrays
  var scores = new Array();
  var residueScore = 3;
  var gapPenalty = -1;
  var misMatchPenalty = -1;
  for(var i = 0; i < seq1.length; ++i) {
    scores[i] = new Array;
    for(var j = 0; j < seq2.length; ++j) {
      scores[i][j] = new Array;
    }
  }

  //set scores
  for(var i = 0; i < seq1.length; ++i) {
    for(var j = 0; j < seq2.length; ++j) {
      for(var k = 0; k < seq3.length; ++k) {
        //Find possible precursor states
        var oldDiag = ( i > 0 && j > 0 && k > 0 ) ? scores[i-1][j-1][k-1].val : 0, //No gap
            oldGap12 = ( k > 0 ) ? scores[i][j][k-1].val : 0, //There is a gap in seq2&3
            oldGap13 = ( j > 0 ) ? scores[i][j-1][k].val : 0, //There is a gap in seq1&3
            oldGap23 = ( i > 0 ) ? scores[i-1][j][k].val : 0, //There is a gap in seq2&1
            oldGap3 = ( i > 0 && j > 0 ) ? scores[i-1][j-1][k].val : 0, //There is a gap in seq3
            oldGap2 = ( i > 0 && k > 0 ) ? scores[i-1][j][k-1].val : 0, //There is a gap in seq2
            oldGap1 = ( j > 0 && k > 0 ) ? scores[i][j-1][k-1].val : 0; //There is a gap in seq1

        //calculate possible states
        var _12 = (seq1[i] == seq2[j]) ? residueScore : misMatchPenalty,
          _13 = (seq1[i] == seq3[k]) ? residueScore : misMatchPenalty,
          _23 = (seq2[j] == seq3[k]) ? residueScore : misMatchPenalty;

        matches_1 = (_12 > 0 || _13 > 0) ? "y" : "n",
          matches_2 = (_12 > 0 || _23 > 0) ? "y" : "n",
          matches_3 = (_13 > 0 || _23 > 0) ? "y" : "n";

        var diag = oldDiag + _12 + _13 + _23;
          gap1 = oldGap1 + gapPenalty + _23,
          gap2 = oldGap2 + gapPenalty + _13,
          gap3 = oldGap3 + gapPenalty + _12,
          gap12 = oldGap12 + 2*gapPenalty,
          gap23 = oldGap23 + 2*gapPenalty,
          gap13 = oldGap13 + 2*gapPenalty;

        //find the best state
        var max = Math.max(diag, gap1, gap2, gap3, gap12, gap23, gap13);

        if( diag >= max ) {
          scores[i][j][k] ={val: diag, point:"diag", matches: [matches_1, matches_2, matches_3]};
        } else if ( gap1 >= max ) {
          scores[i][j][k] ={val: gap1, point:"gap1", matches: ["g", matches_2, matches_3]};
        } else if ( gap2 >= max ) {
          scores[i][j][k] ={val: gap2, point:"gap2", matches: [matches_1, "g", matches_3]};
        } else if( gap3 >= max ) {
          scores[i][j][k] ={val: gap3, point:"gap3", matches: [matches_1, matches_2, "g"]};
        } else if( gap13 >= max ) {
          scores[i][j][k] ={val: gap13, point:"gap13", matches: ["g", "g", "g"]};
        } else if( gap23 >= max ) {
          scores[i][j][k] ={val: gap23, point:"gap23", matches: ["g", "g", "g"]};
        } else if( gap12 >= max ) {
          scores[i][j][k] ={val: gap12, point:"gap12", matches: ["g", "g", "g"]};
        }
      }
    }
  }
  
  //Follow array back
  var i = seq1.length-1;
  var j = seq2.length-1;
  var k = seq3.length-1;
  while(i >= 0 && j >= 0 && k >= 0) {
    var cur = scores[i][j][k];
    aligns_1 = cur.matches[0] + aligns_1;
    aligns_2 = cur.matches[1] + aligns_2;
    aligns_3 = cur.matches[2] + aligns_3;
    switch(cur.point) {
      case "diag":
        s1_a = seq1[i] + s1_a;
        s2_a = seq2[j] + s2_a;
        s3_a = seq3[k] + s3_a;
        --i;
        --j;
        --k;
        break;
      case "gap3":
        s1_a = seq1[i] + s1_a;
        s2_a = seq2[j] + s2_a;
        s3_a = "-" + s3_a;
        --i;
        --j;
        break;
      case "gap2":
        s1_a = seq1[i] + s1_a;
        s2_a = "-"  + s2_a;
        s3_a = seq3[k] + s3_a;
        --i;
        --k;
        break;
      case "gap1":
        s1_a = "-" + s1_a;
        s2_a = seq2[j] + s2_a;
        s3_a = seq3[k] + s3_a;
        --j;
        --k;
        break;
      case "gap13":
        s1_a = "-" + s1_a;
        s2_a = seq2[j] + s2_a;
        s3_a = "-" + s3_a;
        --j;
        break;
      case "gap23":
        s1_a = seq1[i] + s1_a;
        s2_a = "-" + s2_a;
        s3_a = "-" + s3_a;
        --i;
        break;
      case "gap12":
        s1_a = "-" + s1_a;
        s2_a = "-" + s2_a;
        s3_a = seq3[k] + s3_a;
        --k;
        break;
    }
  }

  //Deal with early ends
  while(i >= 0 || j >= 0 || k >= 0) {
    if(i >= 0) {
      s1_a = seq1[i] + s1_a;
      aligns_1 = (seq1[i] == seq2[j] || seq1[i] == seq3[k]) ? "y" : "n" + aligns_1;
      --i;
    } else {
      s1_a = "-" + s1_a;
      aligns_1 = "g" + aligns_1;
    }
    if(j >= 0) {
      s2_a = seq2[j] + s2_a;
      aligns_2 = (seq1[i] == seq2[j] || seq2[j] == seq3[k]) ? "y" : "n" + aligns_2;
      --j;
    } else {
      s2_a = "-" + s2_a;
      aligns_2 = "g" + aligns_2;
    }
    if(k >= 0) {
      aligns_3 = (seq3[k] == seq2[j] || seq1[i] == seq3[k]) ? "y" : "n" + aligns_3;
      s3_a = seq3[k] + s3_a;
      --k;
    } else {
      s3_a = "-" + s3_a;
      aligns_3 = "g" + aligns_3;
    }
  }
}