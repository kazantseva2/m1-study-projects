package org.vmk.dep508.io.rest;

import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;

public class ResponseBuilder {

    private Map<String, String> headers;
    private String body;
    private HttpStatus status;

    public ResponseBuilder() {
        headers = new HashMap<>();
    }

    public ResponseBuilder addHeader(String key, String value) {
        headers.put(key, value);
        return this;
    }

    public ResponseBuilder setBody(String body) {
        this.body = body;
        return this;
    }

    public ResponseBuilder setStatus(HttpStatus status) {
        this.status = status;
        return this;
    }

/*
Пример ответа на запрос
    /*
HTTP/1.1 200
status: 200
content-type: application/json
cache-control: no-cache,no-store,max-age=0,must-revalidate

{"status": "sucess"}


     */
    public void write(OutputStream output) throws IOException {
        if(this.status == null){
            throw new IllegalArgumentException("Please, provide HttpStatus");
        }
        if(this.body == null){
            throw new IllegalArgumentException("Please, provide body");
        }


        // Формируем строку статуса
        StringBuilder response = new StringBuilder();
        response.append("HTTP/1.1 ").append(this.status.getCode()).append("\n");
        response.append("status: ").append(this.status.getCode()).append("\n");

        // Добавляем заголовки
        for (Map.Entry<String, String> header : headers.entrySet()) {
            response.append(header.getKey()).append(": ").append(header.getValue()).append("\n");
        }

        // Пустая строка для разделения заголовков и тела
        response.append("\n");

        // Добавляем тело
        response.append(this.body);

        // Записываем весь ответ в поток
        output.write(response.toString().getBytes());
        output.flush();
    }

    public static void write404(OutputStream output) throws IOException {
        ResponseBuilder builder = new ResponseBuilder();
        builder
                .setStatus(HttpStatus.NOT_FOUND)
                .setBody("<h1>Page not found</h1>")
                .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                .addHeader("content-type", "text/html; charset=UTF-8")
                .write(output);
    }

    public static void writeError(OutputStream output, Exception exception) throws IOException {
        ResponseBuilder builder = new ResponseBuilder();
        builder
                .setStatus(HttpStatus.SERVER_ERROR)
                .setBody("<h1>Server error</h1><br /><pre>"+exception.getMessage()+"</pre>")
                .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                .addHeader("content-type", "text/html; charset=UTF-8")
                .write(output);
    }



    public static void writeSuceess(OutputStream output) throws IOException {
        ResponseBuilder builder = new ResponseBuilder();
        builder
                .setStatus(HttpStatus.OK)
                .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                .addHeader("content-type", "application/json")
                .setBody("{\"status\": \"sucess\"}")
                .write(output);
    }

    public static void writeFailure(OutputStream output) throws IOException {
        ResponseBuilder builder = new ResponseBuilder();
        builder
                .setStatus(HttpStatus.OK)
                .addHeader("cache-control", "no-cache,no-store,max-age=0,must-revalidate")
                .addHeader("content-type", "application/json")
                .setBody("{\"status\": \"failed\"}")
                .write(output);
    }
}
