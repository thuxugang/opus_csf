#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
sys.path.append("..")

import pymongo

class PyConnect(object):
    
    def __init__(self, host, port):
        try:
            self.conn = pymongo.MongoClient(host, port)
        except:
            print ('connect to %s:%s fail' %(host, port))
            exit(0)

    def __del__(self):
        self.conn.close()

    def use(self, dbname):
        self.db = self.conn[dbname]

    def setCollection(self, collection):
        if not self.db:
            print ('don\'t assign database')
            exit(0)
        else:
            self.coll = self.db[collection]

    def find(self, query = {}):
        if type(query) is not dict:
            print ('the type of query isn\'t dict')
            exit(0)
        try:
            if not self.coll:
                print ('don\'t assign collection')
            else:
                result = self.coll.find_one(query)
        except NameError:
            print ('some fields name are wrong in ',query)
            exit(0)
        return result

    def insert(self, data):
        if type(data) is not dict:
            print ('the type of insert data isn\'t dict')
            exit(0)
        self.coll.insert(data)

    def remove(self, data):
        if type(data) is not dict:
            print ('the type of remove data isn\'t dict')
            exit(0)
        self.coll.remove(data)

    def update(self, data, setdata):
        if type(data) is not dict or type(setdata) is not dict:
            print ('the type of update and data isn\'t dict')
            exit(0)
        self.coll.update(data,{'$set':setdata})

def getConnect(db_name,collection_name):
    connect = PyConnect('localhost', 27017)
    connect.use(db_name)
    connect.setCollection(collection_name)
    return connect

