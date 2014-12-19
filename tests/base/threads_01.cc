// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <base/thread_management.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>
template <int> struct X {};
struct U {
  virtual ~U () {}
  X<0> foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_0_const ()const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_0_const () const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_0 () {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_1_const (X<1>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_1_const (X<1>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_1_const (X<1>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_1_const (X<1>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_1_const (X<1>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_1_const (X<1>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_1_const (X<1>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_1_const (X<1>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_1_const (X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_1_const (X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_1_const (const X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_1_const (const X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_1_const (X<1>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_1_const (X<1>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_1_const (X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_1_const (X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_1_const (const X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_1_const (const X<1>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_1 (X<1>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_1 (X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_1 (const X<1>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_2_const (X<1>,X<2>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_2_const (X<1>&,X<2>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_2_const (X<1>&,X<2>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_2_const (X<1>,X<2>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_2_const (X<1>&,X<2>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_2_const (X<1>&,X<2>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_2_const (X<1>,X<2>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_2_const (X<1>,X<2>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_2_const (X<1>&,X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_2_const (X<1>&,X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_2_const (const X<1>&,const X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_2_const (const X<1>&,const X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_2_const (X<1>,X<2>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_2_const (X<1>,X<2>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_2_const (X<1>&,X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_2_const (X<1>&,X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_2_const (const X<1>&,const X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_2_const (const X<1>&,const X<2>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_2 (X<1>,X<2>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_2 (X<1>&,X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_2 (const X<1>&,const X<2>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_3_const (X<1>,X<2>,X<3>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_3_const (X<1>&,X<2>&,X<3>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_3_const (X<1>&,X<2>&,X<3>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_3_const (X<1>,X<2>,X<3>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_3_const (X<1>&,X<2>&,X<3>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_3_const (X<1>&,X<2>&,X<3>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_3_const (X<1>,X<2>,X<3>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_3_const (X<1>,X<2>,X<3>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_3_const (X<1>&,X<2>&,X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_3_const (X<1>&,X<2>&,X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_3_const (const X<1>&,const X<2>&,const X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_3_const (const X<1>&,const X<2>&,const X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_3_const (X<1>,X<2>,X<3>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_3_const (X<1>,X<2>,X<3>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_3_const (X<1>&,X<2>&,X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_3_const (X<1>&,X<2>&,X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_3_const (const X<1>&,const X<2>&,const X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_3_const (const X<1>&,const X<2>&,const X<3>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_3 (X<1>,X<2>,X<3>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_3 (X<1>&,X<2>&,X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_3 (const X<1>&,const X<2>&,const X<3>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_4_const (X<1>,X<2>,X<3>,X<4>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_4_const (X<1>,X<2>,X<3>,X<4>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_4_const (X<1>,X<2>,X<3>,X<4>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_4_const (X<1>,X<2>,X<3>,X<4>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_4_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_4_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_4_const (X<1>,X<2>,X<3>,X<4>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_4_const (X<1>,X<2>,X<3>,X<4>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_4_const (X<1>&,X<2>&,X<3>&,X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_4_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_4_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_4 (X<1>,X<2>,X<3>,X<4>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_4 (X<1>&,X<2>&,X<3>&,X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_4 (const X<1>&,const X<2>&,const X<3>&,const X<4>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_5_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_5_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_5_const (X<1>,X<2>,X<3>,X<4>,X<5>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_5_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_5_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_5_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_5 (X<1>,X<2>,X<3>,X<4>,X<5>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_5 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_5 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_6_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_6_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_6_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_6_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_6_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_6_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_6 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_6 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_6 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_7_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_7_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_7_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_7_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_7_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_7_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_7 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_7 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_7 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  const X<0> & const_ref_foo_const_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual const X<0> & virtual_const_ref_foo_const_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&)const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> static_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static X<0> & static_ref_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  static const X<0> & static_const_ref_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> foo_const_ref_8_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  X<0> & ref_foo_const_ref_8_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_8_const (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_8_const (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_8_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_8_const (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) const {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_8 (X<1>,X<2>,X<3>,X<4>,X<5>,X<6>,X<7>,X<8>) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_ref_8 (X<1>&,X<2>&,X<3>&,X<4>&,X<5>&,X<6>&,X<7>&,X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> virtual_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
  virtual X<0> & virtual_ref_foo_const_ref_8 (const X<1>&,const X<2>&,const X<3>&,const X<4>&,const X<5>&,const X<6>&,const X<7>&,const X<8>&) {
    deallog << __PRETTY_FUNCTION__ << std::endl;
    static X<0> x; return x;
  }
};
int main () {
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  using namespace Threads;
  ThreadGroup<X<0> > tg;
  ThreadGroup<X<0>&> tgr;
  ThreadGroup<const X<0>&> tgcr;
  U u;
X<1> x1;
X<2> x2;
X<3> x3;
X<4> x4;
X<5> x5;
X<6> x6;
X<7> x7;
X<8> x8;
    tgr += spawn (u, &U::ref_foo_0) ();
    tgr += spawn (u, &U::ref_foo_0_const) ();
    tgr += spawn (u, &U::ref_foo_const_ref_0) ();
    tgr += spawn (u, &U::ref_foo_const_ref_0_const) ();
    tgr += spawn (u, &U::ref_foo_ref_0) ();
    tgr += spawn (u, &U::ref_foo_ref_0_const) ();
    tgcr += spawn (u, &U::const_ref_foo_0) ();
    tgcr += spawn (u, &U::const_ref_foo_0_const) ();
    tgcr += spawn (u, &U::const_ref_foo_const_ref_0) ();
    tgcr += spawn (u, &U::const_ref_foo_const_ref_0_const) ();
    tgcr += spawn (u, &U::const_ref_foo_ref_0) ();
    tgcr += spawn (u, &U::const_ref_foo_ref_0_const) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_0) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_0_const) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_0) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_0_const) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_0) ();
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_0_const) ();
    tg += spawn (u, &U::foo_0) ();
    tg += spawn (u, &U::foo_0_const) ();
    tg += spawn (u, &U::foo_const_ref_0) ();
    tg += spawn (u, &U::foo_const_ref_0_const) ();
    tg += spawn (u, &U::foo_ref_0) ();
    tg += spawn (u, &U::foo_ref_0_const) ();
    tgr += spawn (u, &U::virtual_ref_foo_0) ();
    tgr += spawn (u, &U::virtual_ref_foo_0_const) ();
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_0) ();
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_0_const) ();
    tgr += spawn (u, &U::virtual_ref_foo_ref_0) ();
    tgr += spawn (u, &U::virtual_ref_foo_ref_0_const) ();
    tg += spawn (u, &U::virtual_foo_0) ();
    tg += spawn (u, &U::virtual_foo_0_const) ();
    tg += spawn (u, &U::virtual_foo_const_ref_0) ();
    tg += spawn (u, &U::virtual_foo_const_ref_0_const) ();
    tg += spawn (u, &U::virtual_foo_ref_0) ();
    tg += spawn (u, &U::virtual_foo_ref_0_const) ();

    tgr += spawn (&U::static_ref_foo_0) ();
    tgr += spawn (&U::static_ref_foo_const_ref_0) ();
    tgr += spawn (&U::static_ref_foo_ref_0) ();
    tgcr += spawn (&U::static_const_ref_foo_0) ();
    tgcr += spawn (&U::static_const_ref_foo_const_ref_0) ();
    tgcr += spawn (&U::static_const_ref_foo_ref_0) ();
    tg += spawn (&U::static_foo_0) ();
    tg += spawn (&U::static_foo_const_ref_0) ();
    tg += spawn (&U::static_foo_ref_0) ();
    tgr += spawn (u, &U::ref_foo_1) (x1);
    tgr += spawn (u, &U::ref_foo_1_const) (x1);
    tgr += spawn (u, &U::ref_foo_const_ref_1) (x1);
    tgr += spawn (u, &U::ref_foo_const_ref_1_const) (x1);
    tgr += spawn (u, &U::ref_foo_ref_1) (x1);
    tgr += spawn (u, &U::ref_foo_ref_1_const) (x1);
    tgcr += spawn (u, &U::const_ref_foo_1) (x1);
    tgcr += spawn (u, &U::const_ref_foo_1_const) (x1);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_1) (x1);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_1_const) (x1);
    tgcr += spawn (u, &U::const_ref_foo_ref_1) (x1);
    tgcr += spawn (u, &U::const_ref_foo_ref_1_const) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_1) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_1_const) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_1) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_1_const) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_1) (x1);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_1_const) (x1);
    tg += spawn (u, &U::foo_1) (x1);
    tg += spawn (u, &U::foo_1_const) (x1);
    tg += spawn (u, &U::foo_const_ref_1) (x1);
    tg += spawn (u, &U::foo_const_ref_1_const) (x1);
    tg += spawn (u, &U::foo_ref_1) (x1);
    tg += spawn (u, &U::foo_ref_1_const) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_1) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_1_const) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_1) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_1_const) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_ref_1) (x1);
    tgr += spawn (u, &U::virtual_ref_foo_ref_1_const) (x1);
    tg += spawn (u, &U::virtual_foo_1) (x1);
    tg += spawn (u, &U::virtual_foo_1_const) (x1);
    tg += spawn (u, &U::virtual_foo_const_ref_1) (x1);
    tg += spawn (u, &U::virtual_foo_const_ref_1_const) (x1);
    tg += spawn (u, &U::virtual_foo_ref_1) (x1);
    tg += spawn (u, &U::virtual_foo_ref_1_const) (x1);

    tgr += spawn (&U::static_ref_foo_1) (x1);
    tgr += spawn (&U::static_ref_foo_const_ref_1) (x1);
    tgr += spawn (&U::static_ref_foo_ref_1) (x1);
    tgcr += spawn (&U::static_const_ref_foo_1) (x1);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_1) (x1);
    tgcr += spawn (&U::static_const_ref_foo_ref_1) (x1);
    tg += spawn (&U::static_foo_1) (x1);
    tg += spawn (&U::static_foo_const_ref_1) (x1);
    tg += spawn (&U::static_foo_ref_1) (x1);
    tgr += spawn (u, &U::ref_foo_2) (x1,x2);
    tgr += spawn (u, &U::ref_foo_2_const) (x1,x2);
    tgr += spawn (u, &U::ref_foo_const_ref_2) (x1,x2);
    tgr += spawn (u, &U::ref_foo_const_ref_2_const) (x1,x2);
    tgr += spawn (u, &U::ref_foo_ref_2) (x1,x2);
    tgr += spawn (u, &U::ref_foo_ref_2_const) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_2) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_2_const) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_2) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_2_const) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_ref_2) (x1,x2);
    tgcr += spawn (u, &U::const_ref_foo_ref_2_const) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_2) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_2_const) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_2) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_2_const) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_2) (x1,x2);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_2_const) (x1,x2);
    tg += spawn (u, &U::foo_2) (x1,x2);
    tg += spawn (u, &U::foo_2_const) (x1,x2);
    tg += spawn (u, &U::foo_const_ref_2) (x1,x2);
    tg += spawn (u, &U::foo_const_ref_2_const) (x1,x2);
    tg += spawn (u, &U::foo_ref_2) (x1,x2);
    tg += spawn (u, &U::foo_ref_2_const) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_2) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_2_const) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_2) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_2_const) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_ref_2) (x1,x2);
    tgr += spawn (u, &U::virtual_ref_foo_ref_2_const) (x1,x2);
    tg += spawn (u, &U::virtual_foo_2) (x1,x2);
    tg += spawn (u, &U::virtual_foo_2_const) (x1,x2);
    tg += spawn (u, &U::virtual_foo_const_ref_2) (x1,x2);
    tg += spawn (u, &U::virtual_foo_const_ref_2_const) (x1,x2);
    tg += spawn (u, &U::virtual_foo_ref_2) (x1,x2);
    tg += spawn (u, &U::virtual_foo_ref_2_const) (x1,x2);

    tgr += spawn (&U::static_ref_foo_2) (x1,x2);
    tgr += spawn (&U::static_ref_foo_const_ref_2) (x1,x2);
    tgr += spawn (&U::static_ref_foo_ref_2) (x1,x2);
    tgcr += spawn (&U::static_const_ref_foo_2) (x1,x2);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_2) (x1,x2);
    tgcr += spawn (&U::static_const_ref_foo_ref_2) (x1,x2);
    tg += spawn (&U::static_foo_2) (x1,x2);
    tg += spawn (&U::static_foo_const_ref_2) (x1,x2);
    tg += spawn (&U::static_foo_ref_2) (x1,x2);
    tgr += spawn (u, &U::ref_foo_3) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_3_const) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_const_ref_3) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_const_ref_3_const) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_ref_3) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_ref_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_3) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_3) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_ref_3) (x1,x2,x3);
    tgcr += spawn (u, &U::const_ref_foo_ref_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_3) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_3) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_3_const) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_3) (x1,x2,x3);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_3_const) (x1,x2,x3);
    tg += spawn (u, &U::foo_3) (x1,x2,x3);
    tg += spawn (u, &U::foo_3_const) (x1,x2,x3);
    tg += spawn (u, &U::foo_const_ref_3) (x1,x2,x3);
    tg += spawn (u, &U::foo_const_ref_3_const) (x1,x2,x3);
    tg += spawn (u, &U::foo_ref_3) (x1,x2,x3);
    tg += spawn (u, &U::foo_ref_3_const) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_3) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_3_const) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_3) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_3_const) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_ref_3) (x1,x2,x3);
    tgr += spawn (u, &U::virtual_ref_foo_ref_3_const) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_3) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_3_const) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_const_ref_3) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_const_ref_3_const) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_ref_3) (x1,x2,x3);
    tg += spawn (u, &U::virtual_foo_ref_3_const) (x1,x2,x3);

    tgr += spawn (&U::static_ref_foo_3) (x1,x2,x3);
    tgr += spawn (&U::static_ref_foo_const_ref_3) (x1,x2,x3);
    tgr += spawn (&U::static_ref_foo_ref_3) (x1,x2,x3);
    tgcr += spawn (&U::static_const_ref_foo_3) (x1,x2,x3);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_3) (x1,x2,x3);
    tgcr += spawn (&U::static_const_ref_foo_ref_3) (x1,x2,x3);
    tg += spawn (&U::static_foo_3) (x1,x2,x3);
    tg += spawn (&U::static_foo_const_ref_3) (x1,x2,x3);
    tg += spawn (&U::static_foo_ref_3) (x1,x2,x3);
    tgr += spawn (u, &U::ref_foo_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_4_const) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_const_ref_4_const) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_ref_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_ref_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::const_ref_foo_ref_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_4_const) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_const_ref_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_const_ref_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_ref_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::foo_ref_4_const) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_4_const) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_4_const) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_ref_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::virtual_ref_foo_ref_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_const_ref_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_const_ref_4_const) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_ref_4) (x1,x2,x3,x4);
    tg += spawn (u, &U::virtual_foo_ref_4_const) (x1,x2,x3,x4);

    tgr += spawn (&U::static_ref_foo_4) (x1,x2,x3,x4);
    tgr += spawn (&U::static_ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgr += spawn (&U::static_ref_foo_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (&U::static_const_ref_foo_4) (x1,x2,x3,x4);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_4) (x1,x2,x3,x4);
    tgcr += spawn (&U::static_const_ref_foo_ref_4) (x1,x2,x3,x4);
    tg += spawn (&U::static_foo_4) (x1,x2,x3,x4);
    tg += spawn (&U::static_foo_const_ref_4) (x1,x2,x3,x4);
    tg += spawn (&U::static_foo_ref_4) (x1,x2,x3,x4);
    tgr += spawn (u, &U::ref_foo_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_5_const) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_ref_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::const_ref_foo_ref_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_const_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::foo_ref_5_const) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_5_const) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::virtual_ref_foo_ref_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_const_ref_5_const) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (u, &U::virtual_foo_ref_5_const) (x1,x2,x3,x4,x5);

    tgr += spawn (&U::static_ref_foo_5) (x1,x2,x3,x4,x5);
    tgr += spawn (&U::static_ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (&U::static_ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (&U::static_const_ref_foo_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tgcr += spawn (&U::static_const_ref_foo_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (&U::static_foo_5) (x1,x2,x3,x4,x5);
    tg += spawn (&U::static_foo_const_ref_5) (x1,x2,x3,x4,x5);
    tg += spawn (&U::static_foo_ref_5) (x1,x2,x3,x4,x5);
    tgr += spawn (u, &U::ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_6_const) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::const_ref_foo_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::foo_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_6_const) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::virtual_ref_foo_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_const_ref_6_const) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (u, &U::virtual_foo_ref_6_const) (x1,x2,x3,x4,x5,x6);

    tgr += spawn (&U::static_ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (&U::static_ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (&U::static_ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (&U::static_const_ref_foo_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tgcr += spawn (&U::static_const_ref_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (&U::static_foo_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (&U::static_foo_const_ref_6) (x1,x2,x3,x4,x5,x6);
    tg += spawn (&U::static_foo_ref_6) (x1,x2,x3,x4,x5,x6);
    tgr += spawn (u, &U::ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::const_ref_foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::virtual_ref_foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_const_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (u, &U::virtual_foo_ref_7_const) (x1,x2,x3,x4,x5,x6,x7);

    tgr += spawn (&U::static_ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (&U::static_ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (&U::static_ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (&U::static_const_ref_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgcr += spawn (&U::static_const_ref_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (&U::static_foo_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (&U::static_foo_const_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tg += spawn (&U::static_foo_ref_7) (x1,x2,x3,x4,x5,x6,x7);
    tgr += spawn (u, &U::ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::ref_foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::ref_foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::ref_foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::const_ref_foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (u, &U::virtual_ref_foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_const_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (u, &U::virtual_foo_ref_8_const) (x1,x2,x3,x4,x5,x6,x7,x8);

    tgr += spawn (&U::static_ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (&U::static_ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgr += spawn (&U::static_ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (&U::static_const_ref_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (&U::static_const_ref_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tgcr += spawn (&U::static_const_ref_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (&U::static_foo_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (&U::static_foo_const_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
    tg += spawn (&U::static_foo_ref_8) (x1,x2,x3,x4,x5,x6,x7,x8);
  tg.join_all();
  tgr.join_all();
  tgcr.join_all();

  deallog.detach ();
  logfile.close ();
  unify_pretty_function ("output");
  sort_file_contents ("output");
}
