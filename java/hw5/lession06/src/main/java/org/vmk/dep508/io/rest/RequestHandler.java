package org.vmk.dep508.io.rest;

import org.apache.log4j.Logger;

import java.io.*;
import java.net.Socket;
import java.util.Arrays;

public class RequestHandler implements Runnable {

    private Logger log = Logger.getLogger(getClass());
    private Socket socket;

    public RequestHandler(Socket socket) {
        log.debug("New request handler");
        this.socket = socket;
    }

    @Override
    public void run() {
        try {
            SerializationHelper<Book> bookSerializationHelper = new SerializationHelper<>(Book.class);

            Request request = Request.parse(socket.getInputStream());
            ResponseBuilder builder = new ResponseBuilder();

            if(request.getPath().equals("/books")) {

                if(request.getMethod() != HttpMethod.GET) {
                    ResponseBuilder.writeError(socket.getOutputStream(), new Exception("POST method not supported for /books"));
                }
                log.info("Get book list");

                File f = new File("books");
                String list = Arrays
                        .stream(f.listFiles())
                        .map(s -> "/books/"+s.getName())
                        .map(s -> "\""+s+"\"")
                        .reduce((s1, s2) -> s1+", "+s2).get();

                builder
                        .setStatus(HttpStatus.OK)
                        .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                        .addHeader("content-type", "application/json")
                        //.addHeader("content-type", "text/html; charset=UTF-8")
                        .setBody("["+list+"]")
                        .write(socket.getOutputStream());

            } else if(request.getPath().startsWith("/books")) {
                String relativePath = request.getPath().substring(1);

                if(request.getMethod() == HttpMethod.GET) {
                    log.info("Get book");

                    Book result = bookSerializationHelper.loadFromFile(relativePath);
                    if(result != null) {
                        builder
                                .setStatus(HttpStatus.OK)
                                .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                                .addHeader("content-type", "application/json")
                                .setBody(bookSerializationHelper.convertToJsonString(result))
                                .write(socket.getOutputStream());

                    } else {
                        ResponseBuilder.write404(socket.getOutputStream());
                    }

                } else if(request.getMethod() == HttpMethod.POST ) {
                    log.info("Save book");

                    Book toSave = bookSerializationHelper.parseJson(request.getBody());
                    boolean result = bookSerializationHelper.saveToFile(relativePath, toSave);
                    if(result) {
                        ResponseBuilder.writeSuceess(socket.getOutputStream());

                    } else {
                        ResponseBuilder.writeFailure(socket.getOutputStream());
                    }
                } else if(request.getMethod() == HttpMethod.DELETE ) {
                    log.info("Remove book");

                    File bookFile = new File(relativePath);
                    if(bookFile.delete()) {
                        ResponseBuilder.writeSuceess(socket.getOutputStream());
                    } else {
                        ResponseBuilder.writeFailure(socket.getOutputStream());
                    }
                }
            }
            else {
                log.info("Processing 404 page");
                ResponseBuilder.write404(socket.getOutputStream());
            }
            log.info("Response sent");
            socket.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
