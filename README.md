# ska-based LPPM

## Introduction

This repository contains the code and simulated data developed for my Master’s Thesis, "**Protecting sensitive location information for mobile fitness applications with spatial k-anonymity**". The main focus of the thesis is on the development and evaluation of a ska-based geomasking method to protect the privacy of trajectories.

### Contents
- A Python script implementing the developed geomasking method
- A Python script deleting BOM encoding in a GPX file
- Simulated trajectories (GPX format) in four districts in Vienna used for testing and demonstrating the method.

## Usage

To run the script for the ska-based geomasking method you need the following:
- Install the necessary python libraries
- Path to the input GPX file
- Path to the output GPX file
- Define the ska level (k-min and k-max)
- Openrouteservices API key
