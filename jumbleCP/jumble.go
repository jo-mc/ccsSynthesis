// Copyright ©2020 J McConnell	. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Generate synthetic CCS reads from a alignment.
//

package main

import (
	"bufio"
	"fmt"
	"log"
	"math/rand"
	"os"
	"time"
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
	if len(os.Args) <= 1 {
		fmt.Printf("USAGE : %s <target_filename> \n", os.Args[0])
		os.Exit(0)
	}

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
	}

	var fasta []string
	var aLine int

	fdata, e := Readln(r) // my fasta has 3 non fasta format lines at head (also a couple at tail)
	aLine++
	fdata, e = Readln(r)
	aLine++
	fdata, e = Readln(r)
	aLine++
	for { // j := 1; j <= ; j++                ENDLESS LOOP
		fdata, e = Readln(r)
		aLine++
		if e == nil {
			if fdata[0:1] == ">" {
				fasta = append(fasta, fdata)
				fdata, ee := Readln(r) // now get read data
				aLine++
				if ee == nil {
					fasta = append(fasta, fdata)
				} else {
					// not valid fasta ?
					fmt.Printf("error reading fasta format at line %d : first 20 chars of line %s", aLine, fdata[0:20])
					break
				}
			} else {
				// not valid fasta ?
				fmt.Printf("error reading fasta format at line %d : first 20 chars of line %s", aLine, fdata[0:20])
				break
			}

		} else {
			break
		}

	}
	fmt.Println(len(fasta), cap(fasta))

	rand.Seed(time.Now().UnixNano()) // https://golang.cafe/blog/golang-random-number-generator.html

	if len(fasta) > 3 {
		rand.Shuffle(len(fasta)-1, func(i, j int) { // https://yourbasic.org/golang/shuffle-slice-array/
			if (i%2 == 0) && (j%2 == 0) {
				fasta[i], fasta[j] = fasta[j], fasta[i]
				fasta[i+1], fasta[j+1] = fasta[j+1], fasta[i+1]
			}
			if (i%2 == 1) && (j%2 == 1) {
				fasta[i-1], fasta[j-1] = fasta[j-1], fasta[i-1]
				fasta[i], fasta[j] = fasta[j], fasta[i]
			}
			if (i%2 == 0) && (j%2 == 1) {
				fasta[i], fasta[j-1] = fasta[j-1], fasta[i]
				fasta[i+1], fasta[j] = fasta[j], fasta[i+1]
			}
			if (i%2 == 1) && (j%2 == 0) {
				fasta[i-1], fasta[j] = fasta[j], fasta[i-1]
				fasta[i], fasta[j+1] = fasta[j+1], fasta[i]
			}
		})
	}
	// output to file    https://www.golangprograms.com/write-string-slice-line-by-line-to-a-text-file.html
	outfile := "jumble" + os.Args[1]

	file, err := os.OpenFile(outfile, os.O_CREATE|os.O_WRONLY, 0644)

	if err != nil {
		log.Fatalf("failed creating file: %s, %s", err, outfile)
		os.Exit(0)
	}

	datawriter := bufio.NewWriter(file)

	for _, data := range fasta {
		_, _ = datawriter.WriteString(data + "\n")
	}

	datawriter.Flush()
	file.Close()

}
