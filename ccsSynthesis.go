// Copyright ©2020 J McConnell	. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Generate synthetic CCS reads from a alignment.
// COMMAND:  ./ccsSynthesis /PATH/GRCH38chr8-REF.fa > /PATH/chr8_revcompSim.fa
//  argument = region/contig/chromosome to be split into reads
// read length to generate: line 51 ccs =
// read depth : line 52 divisor of ccs   ( ie ccs /20  +> read depth of twenty )

package main

import (
	"bufio"
	"fmt"
	"math"
	"runtime"
	"strconv"
	"time"

	"log"
	"os"
	//"REGSPLT/rs"
	//	ars "regsplt/rs"
)

// Readln returns a single line (without the ending \n)
// from the input buffered reader.
// An error is returned iff there is an error with the
// buffered reader.
// https://stackoverflow.com/questions/6141604/go-readline-string
func Readln(r *bufio.Reader) (string, error) {
	var (
		isPrefix bool  = true
		err      error = nil
		line, ln []byte
	)
	for isPrefix && err == nil {
		line, isPrefix, err = r.ReadLine()
		ln = append(ln, line...)
	}
	return string(ln), err
}

var errStream = os.Stderr
var errFile string
var err error

func main() {

	outfile, err := os.OpenFile("out.fa", os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	defer outfile.Close()
	w := bufio.NewWriter(outfile)

	file, err := os.OpenFile("logs.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	log.SetOutput(file)

	stTime := time.Now()
	log.Println("Starting time is", stTime)
	log.Println("Generating CCS 15 Kb reads 20x")
	ccs := 15000    //15000   test 50
	gap := ccs / 20 // 15000 / 20  test 5
	nccs := ccs
	ngap := gap
	RevCompAlternate := true //false

	pipeorfileOK() // test we have a pipe in or a file in, will quit if neither or both.
	//  OPEN READER file or PIPE!!
	var r *bufio.Reader

	//var file *os.File
	if len(os.Args) > 1 {
		file, err := os.Open(os.Args[1])
		if err != nil {
			log.Println("error opening file= ", err)
			os.Exit(1)
		}
		defer file.Close()
		r = bufio.NewReader(file)
		log.Println("Input file is ", os.Args[1])
	} else {
		r = bufio.NewReader(os.Stdin)
		log.Println("Getting data from PIPE")
	}

	var readDepth int
	var readSize int
	var noVAR string
	noVAR = ""
	if len(os.Args) > 3 {
		readSize, err = strconv.Atoi(os.Args[3])
		if err == nil {
			ccs = readSize
			nccs = ccs
			gap = ccs / 20
			ngap = gap
		}
	}
	if len(os.Args) > 2 {
		readDepth, err = strconv.Atoi(os.Args[2])
		if err == nil {
			gap = ccs / readDepth
			ngap = gap
		}

	}
	if len(os.Args) > 4 {
		noVAR = os.Args[4]
	}
	var linesRead uint32 = 0
	var regions int = 0
	var index int = 1
	var rID int = 0
	var bp int = 0
	var aRead string
	var rName string
	var fasta []rune
	for { // j := 1; j <= ; j++                ENDLESS LOOP
		fdata, e := Readln(r)
		linesRead++
		if e == nil && len(fdata) > 0 { //  and not blank line?  will exit if have a blank line?  -> loop over till line with data or eof?
			if fdata[0:1] == ">" {
				rName = fdata
				regions++
			} else {
				bp += len(fdata)
				fasta = append(fasta, []rune(fdata)...)
				for len(fasta) >= nccs { // while loop!
					//aRead = string(fasta[:nccs])
					aRead = string(fasta[22 : nccs-22]) // have 22bp gap b.w reads
					rID += 1
					if rID == 1 {
						fmt.Fprintf(w, ">I%09dI%03d:%05d:I%dI 22bp b/w gap:ccs:index: Settings- gap %d target size: %d %s \n", rID, ngap, nccs, index, gap, ccs, rName)
					} else {
						if RevCompAlternate {
							fmt.Fprintf(w, ">I%09dI%03d:%05d:I%dIRevCompI%s\n", rID, ngap, nccs, index, rName)
						} else {
							fmt.Fprintf(w, ">I%09dI%03d:%05d:I%dI%s\n", rID, ngap, nccs, index, rName)
						}
					}
					if RevCompAlternate {
						if rID%2 == 0 {
							fmt.Fprintf(w, "%s\n", aRead)
						} else {
							fmt.Fprintf(w, "%s\n", revComp(aRead))
						}
					} else {
						fmt.Fprintf(w, "%s\n", aRead)
					}
					// variation
					theta := math.Round(math.Mod(float64(rID), 31)) / 10  // 0..3.1
					ngap = gap + int((math.Sin(theta*2))*float64(gap)/4)  //  gap +/- gap/4
					nccs = ccs + int((math.Sin(theta*4))*float64(ccs)/10) //  gap +/- gap/4
					if noVAR == "novar" {
						// no variation in gap or read length
						ngap = gap
						nccs = ccs
					}

					index += ngap
					fasta = fasta[ngap:]
				}
			}
		} else {
			log.Println("end of input: Regions processed: ", regions)
			//fmt.Println("end of input: : ", regions)
			log.Println("reads Written: ", rID)
			log.Println(" base pairs processed: ", bp)    //  index+ccs will be close
			log.Println(" base pairs written: ", ccs*rID) //  index+ccs will be close
			// fmt.Printf(">|%09d|%s\n", rID, rName)
			// fmt.Printf("%s\n", aRead)
			break
		}

	}
	log.Println("Ending time is", time.Now(), " Lines read: ", linesRead)
}

func revComp(aR string) (x string) {
	aRev := []rune(aR)
	// reverse
	for i := len(aRev)/2 - 1; i >= 0; i-- {
		opp := len(aRev) - 1 - i
		aRev[i], aRev[opp] = aRev[opp], aRev[i]
	}
	// complement
	for i := 0; i < len(aRev); i++ {
		aBP := string(aRev[i])
		switch aBP {
		case "A":
			aRev[i] = 'T'
		case "T":
			aRev[i] = 'A'
		case "G":
			aRev[i] = 'C'
		case "C":
			aRev[i] = 'G'
		case "a":
			aRev[i] = 't'
		case "t":
			aRev[i] = 'a'
		case "g":
			aRev[i] = 'c'
		case "c":
			aRev[i] = 'g'
			//default :
		}
	}
	x = string(aRev)
	return
}

func pipeorfileOK() {

	var info os.FileInfo
	var err error
	var inPipe bool = false
	var inFile bool = false

	if runtime.GOOS == "windows" {
		log.Println("- -   -  - > Windows detected. Note: Window pipes not implemented, file argument ok.")
	} else {
		// Do we have piped input?
		info, err = os.Stdin.Stat() // standard input file descriptor
		if err != nil {
			fmt.Println("error reading stdin - exiting")
			panic(err)
		}
		if info.Mode()&os.ModeNamedPipe != 0 { // is data begin piped in?
			// we have a pipe input
			log.Println("we have a pipe input")
			inPipe = true
		}
	}

	// Do we have argument input?
	//var file *os.File
	if len(os.Args) > 1 { // do we have arguments : ie a file to read?
		log.Print(os.Args[1])
		log.Println(" : argument (file) input")
		file, err := os.Open(os.Args[1])
		if err != nil {
			fmt.Println("error opening file= ", err)
			os.Exit(1)
		}
		file.Close()
		inFile = true
	}

	if runtime.GOOS != "windows" {
		// Both pipe and argument? -> EXIT
		if inPipe && inFile {
			fmt.Println("- -   -  - > we have a pipe input and a file input ?? Please use one only, exiting")
			os.Exit(1)
		}
	}

	if (inPipe || inFile) == false {
		// no input
		fmt.Println("- -   ->>->> No input detected ?? exiting. Create reads from a genome (fasta formatted)")
		fmt.Println("Generate synthetic CCS reads from a alignment. Usage:")
		fmt.Printf("\n             ccsSynthesis /PATH/GRCH38chr8-REF.fa > /PATH/chr8_revcompSim.fa\n\n")
		fmt.Println("- -   -  - > Usage: Pipe numbers into program (Linux only): awk '{ print $0 }' datafile.fa I ccsSynthesis")
		fmt.Println("- -   -  - > or use with a file argument (Linux or Windows): ccsSynthesis datafile.fa")
		fmt.Println("- -   -  - > File input must be in fasta format.")
		fmt.Println("- -   -  - > Alter by appending required values, read depth and (optionally) read length to input file ")
		fmt.Println("- -   -  - > ./ccsSynthesis /PATH/GRCH38chr8-REF.fa 10 10000 > /PATH/chr8_revcompSim.fa")
		fmt.Println("- -   -  - > can have new read depth without read length (15000 will be used).")
		fmt.Println("- -   -  - > with 'novar' appended to depth and length supplied values, fixed lengths and gaps will be used (no variation).")
		fmt.Printf("\n- -   -  - > ccsSynthesis datafile.fa 20 10000 novar      With no parameters after file the Default settings are:  15,000 read length and 20x Coverage, variation.\n")
		os.Exit(1)
	}
}
