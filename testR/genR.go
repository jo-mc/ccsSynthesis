// Copyright ©2020 J McConnell	. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Generate synthetic CCS reads from a alignment.

// Test : takes first read/region of supplied fasta file and generates four files:
//  the first read ref file ref.fa, a reversed rev.fa, a complemented comp.fa and a rev complement revcomp.fa
//

package main

import (
	"bufio"
	"fmt"
	"runtime"
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

	file, err := os.OpenFile("logs.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}

	log.SetOutput(file)

	reffile, err := os.OpenFile("ref.fa", os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(reffile)

	revfile, err := os.OpenFile("comp.fa", os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	x := bufio.NewWriter(revfile)

	compfile, err := os.OpenFile("rev.fa", os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	y := bufio.NewWriter(compfile)

	revcompfile, err := os.OpenFile("revcomp.fa", os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	z := bufio.NewWriter(revcompfile)

	stTime := time.Now()
	log.Println("Test generation", stTime)
	log.Println("Generating ref, complement, reverse and reverse complement for supplied file first read.")

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

	var aRead string
	var rName string
	var fasta []rune
	var aRev []rune
	var maxSize int = 51000

	fdata, e := Readln(r) // get read nmae
	if e == nil {
		rName = fdata
	} else {
		fmt.Println("error reading file= ", err)
	}

	fdata, e = Readln(r) // get 1st fasta sequence (check it is long enough for what you want!)
	if e != nil {
		fmt.Println("error reading file= ", err)
	}
	fasta = append(fasta, []rune(fdata)...)
	if len(fasta) > maxSize {
		fasta = fasta[1:maxSize]
	}
	aRead = string(fasta)
	fmt.Fprintf(w, ">REF:%s\n", rName)
	fmt.Fprintf(w, "%s\n", aRead)

	// complement
	//aRev = make([]string, 0,1)
	aRev = []rune(aRead)
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
	fmt.Fprintf(x, ">COMP:%s\n", rName)
	fmt.Fprintf(x, "%s\n", string(aRev))

	// reverse
	//aRev = make([]string, 0,1)
	aRev = []rune(aRead)
	for i := len(aRev)/2 - 1; i >= 0; i-- {
		opp := len(aRev) - 1 - i
		aRev[i], aRev[opp] = aRev[opp], aRev[i]
	}
	fmt.Fprintf(y, ">Rev:%s\n", rName)
	fmt.Fprintf(y, "%s\n", string(aRev))

	// Reverse complement
	//aRev = make([]string, 0,1)
	//aRev = []rune(aRead)
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
	fmt.Fprintf(z, ">RevComp:%s\n", rName)
	fmt.Fprintf(z, "%s\n", string(aRev))

	log.Println("Ending time is", time.Now())
	w.Flush()
	x.Flush()
	y.Flush()
	z.Flush()

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
			//			fmt.Println("we have a pipe input")
			inPipe = true
		}
	}

	// Do we have argument input?
	//var file *os.File
	if len(os.Args) > 1 { // do we have arguments : ie a file to read?
		//		fmt.Print(os.Args[1])
		//		fmt.Println(" : argument (file) input")
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
