// Copyright ©2020 J McConnell	. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Generate synthetic CCS reads from a alignment.
//

package main

import (
	"bufio"
	"fmt"
	"math"
	"runtime"
	"time"

	//"log"
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
	stTime := time.Now()
	fmt.Println("Starting time is", stTime)
	fmt.Println("Generating CCS 15 Kb reads 20x")
	ccs := 15000    //15000   test 50
	gap := ccs / 20 // 15000 / 20  test 5
	nccs := ccs
	ngap := gap

	pipeorfileOK() // test we have a pipe in or a file in, will quit if neither or both.
	//  OPEN READER file or PIPE!!
	var r *bufio.Reader

	//var file *os.File
	if len(os.Args) > 1 {
		file, err := os.Open(os.Args[1])
		if err != nil {
			fmt.Println("error opening file= ", err)
			os.Exit(1)
		}
		defer file.Close()
		r = bufio.NewReader(file)
	} else {
		r = bufio.NewReader(os.Stdin)
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
		if e == nil {
			if fdata[0:1] == ">" {
				rName = fdata
				regions++
			} else {
				bp += len(fdata)
				fasta = append(fasta, []rune(fdata)...)
				for len(fasta) >= nccs { // while loop!
					aRead = string(fasta[:nccs])
					rID += 1
					if rID == 1 {
						fmt.Printf(">|%09d|%03d:%05d:|%d| gap:ccs:index: Settings- gap %d target size: %d %s \n", rID, ngap, nccs, index, gap, ccs, rName)
					} else {
						fmt.Printf(">|%09d|%03d:%05d:|%d|%s\n", rID, ngap, nccs, index, rName)
					}
					//if rID%2 == 0 {
					fmt.Printf("%s\n", aRead)
					//} else {
					//	fmt.Printf("%s\n", revComp(aRead))
					//}
					// variation
					theta := math.Round(math.Mod(float64(rID), 31)) / 10  // 0..3.1
					ngap = gap + int((math.Sin(theta*2))*float64(gap)/4)  //  gap +/- gap/4
					nccs = ccs + int((math.Sin(theta*4))*float64(ccs)/10) //  gap +/- gap/4
					// no variation in gap or read length
					// ngap = gap
					// nccs = ccs

					index += ngap
					fasta = fasta[ngap:]
				}

			}
		} else {
			fmt.Println("end of input: Regions processed: ", regions)
			//fmt.Println("end of input: : ", regions)
			fmt.Println("reads Written: ", rID)
			fmt.Println(" base pairs processed: ", bp)    //  index+ccs will be close
			fmt.Println(" base pairs written: ", ccs*rID) //  index+ccs will be close
			// fmt.Printf(">|%09d|%s\n", rID, rName)
			// fmt.Printf("%s\n", aRead)
			break
		}

	}
	fmt.Println("Ending time is", time.Now(), " Lines read: ", linesRead)
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
		fmt.Println("- -   -  - > Windows detected. Note: Window pipes not implemented, file argument ok.")
	} else {
		// Do we have piped input?
		info, err = os.Stdin.Stat() // standard input file descriptor
		if err != nil {
			fmt.Println("error reading stdin - exiting")
			panic(err)
		}
		if info.Mode()&os.ModeNamedPipe != 0 { // is data begin piped in?
			// we have a pipe input
			fmt.Println("we have a pipe input")
			inPipe = true
		}
	}

	// Do we have argument input?
	//var file *os.File
	if len(os.Args) > 1 { // do we have arguments : ie a file to read?
		fmt.Print(os.Args[1])
		fmt.Println(" : argument (file) input")
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
		fmt.Println("- -   -  - > No input detected ?? exiting")
		fmt.Println("- -   -  - > Usage: Pipe numbers into program (Linux only)")
		fmt.Println("- -   -  - > awk '{ print $3 }' datafile.dat | nebulostat")
		fmt.Println("- -   -  - > or use with a file argument (Linux or Windows)")
		fmt.Println("- -   -  - > nebulostat datafile.dat")
		fmt.Println("- -   -  - > or awk version")
		fmt.Println("- -   -  - > awk -f nebulostat.awk datafile.dat,   ")
		fmt.Println("- -   -  - > or pipe in:")
		fmt.Println("- -   -  - > awk '{ print $3 }' datafile.dat | awk -f nebulostat.awk")
		fmt.Println("- -   -  - > File input must consist of one number per line.")
		os.Exit(1)
	}
}
