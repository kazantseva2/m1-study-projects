package org.vmk.dep508.io.rest;

import org.apache.log4j.Logger;

import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Server {

    private static Logger log = Logger.getLogger(Server.class);

    public static void main(String[] args) throws IOException {

        int port = 8080;

        log.info("Starting server on " + port + " port");
        ServerSocket server = new ServerSocket(port);

        ExecutorService service = Executors.newFixedThreadPool(10);

        log.info("Waiting for client connection");
        while(true) {
            Socket s = server.accept();
            log.info("Client connected to server");

            service.submit(new RequestHandler(s));
        }
    }
}
